#!/usr/bin/env python3

# Todo:
# download galaxy catalog
# map names to HIP 
# nonfixtarget
# query objects from internet database and cache them to a catalog file
# speed up calculations
# load list of Planets
# load list of known stars
# load list of messier objects
# filter stars by mag
# filter messier by: mag, type
# for a given object in stars or messier or planets: produce
#     plot the motion in the sky for a given hour of the night with passing of days
# for a given date and a list of objects: produce
#     filter the list based on observability for that date
#     plot their motion in the sky


import astroplan, astroplan.plots
import astropy.coordinates, astropy.time
import astropy.units as u
import numpy,pandas, pytz,datetime,click
import matplotlib.pyplot as plt
import astroquery.vizier

astropy.coordinates.solar_system_ephemeris.set('jpl')


polse=astroplan.Observer(
    latitude=
        46+(28+9.1/60)/60,
    longitude=
        13+(0+42.6/60)/60,
    elevation=
        750*u.meter,
    timezone=
        pytz.timezone('Europe/Rome'))
        
polse_constraints = [
    astroplan.AltitudeConstraint(min=20*u.deg),
    astroplan.AtNightConstraint.twilight_civil(),
    astroplan.MoonSeparationConstraint(min=20*u.deg),
    astroplan.LocalTimeConstraint(min=datetime.time(20,0),max=datetime.time(1,0))
    ]
    
months_labels=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

class Catalog():    
    '''
    My favorite objects must be loaded into
    a catalog.
    '''
    def __init__(self,filename):
        self.data = pandas.read_csv(filename,index_col='name')
        
catalog = Catalog('catalog.csv')


@click.group()
def cli():
    pass

def target_by_name(name,time):
    if name in astropy.coordinates.solar_system_ephemeris.bodies:
        target=astropy.coordinates.get_body(name,time)
        target.name=name # ugly fix to a problem in astropy
    else:
        target = astroplan.FixedTarget.from_name(name)
    return target


class rise_time():
    '''
    Functional class, that computes for
    a given time point, the time it takes for
    the object to rise, relative to the observer.
    '''
    def __init__(self,observer,name,label):
        self.name = name
        self.observer = observer 
        self.label = label
        if name in astropy.coordinates.solar_system_ephemeris.bodies:
            self.isfixed=False
            self.target = None
        else:
            self.isfixed=True
            self.target = astroplan.FixedTarget.from_name(name)
    def __call__(self,time):
        if self.name=='Sun':
            dt = self.observer.sun_rise_time(
                time,which='next')-time
        else:
            if not self.isfixed:
                self.target = target_by_name(self.name,time)
            dt = self.observer.target_rise_time(
                time,self.target,which='next')-time
        return dt.to(u.hour).value

class set_time():
    '''
    Functional class, that computes for
    a given time point, the time it takes (in hours) for the
    object to set, relative to the observer.
    '''
    def __init__(self,observer,name,label):
        self.name = name
        self.observer = observer 
        self.label = label
        if name in astropy.coordinates.solar_system_ephemeris.bodies:
            self.isfixed=False
            self.target = None
        else:
            self.isfixed=True
            self.target = astroplan.FixedTarget.from_name(name)
    def __call__(self,time):
        if self.name=='Sun':
            dt = self.observer.sun_set_time(
                time,which='next')-time
        else:
            if not self.isfixed:
                self.target = target_by_name(self.name,time)
            dt = self.observer.target_set_time(
                time,self.target,which='next')-time
        return dt.to(u.hour).value

class observability():
    '''
    Functional class, that computes for
    a given time point, the observability of the
    object within the following 24 hours, 
    relative to the observer and
    defined constraints.
    '''
    def __init__(self,observer,name,constraints,label):
        self.name = name
        self.observer = observer 
        self.label = label
        self.constraints = constraints
        if name in astropy.coordinates.solar_system_ephemeris.bodies:
            self.isfixed=False
            self.target = None
        else:
            self.isfixed=True
            self.target = astroplan.FixedTarget.from_name(name)
    def __call__(self,time,time_grid_resolution=10*u.minute):
        if self.name=='Sun':
            return numpy.nan
        t1 = time
        t2 = time + 1*u.day
        if not self.isfixed:
            self.target = target_by_name(self.name,time)
        ot = astroplan.observability_table(
            self.constraints,self.observer,[self.target],
                time_range = astropy.time.Time([t1,t2]),
                time_grid_resolution = time_grid_resolution)
        otime = ot[0]['fraction of time observable']*u.day
        return otime.to(u.hour).value
        
class ephemeris():
    def __init__(self,observer,name,fields):
        self.name=name
        self.observer = observer
        self.fields = fields
        self.data = None
        try:
            self.data = pandas.read_csv(self.name + '.eph',parse_dates=['date'],index_col='date')
        except:
            self.data = pandas.DataFrame( 
                columns = ['date'] + [f.label for f in self.fields]  ).set_index('date')
            self.compute(pandas.date_range(start='2020-01-01',end='2020-12-31'))
  
    def compute(self,date_range):
        data_dict = { 'date': list()  }
        for f in self.fields:
            data_dict[f.label]=list()
        with click.progressbar(date_range,label="Computing %s ephemeris" % self.name) as bar:
            for d in bar:
                if d in self.data.index:
                    continue
                data_dict['date'].append(d)
                tz = self.observer.timezone
                dt = datetime.datetime.combine(d.date(),datetime.time(12,0))
                T = astropy.time.Time(tz.localize(dt))
                for f in self.fields:
                    data_dict[f.label].append( f(T))
        df = pandas.DataFrame(data_dict).set_index('date')
        self.data = self.data.append(df)
        
    def __del__(self):
        self.data.to_csv( self.name+".eph")
    

def get_observability(name,dates):
    eph = ephemeris(
        polse,
        name,
        [   
            rise_time(polse,name,'rise'),
            set_time(polse,name,'set'),
            observability(polse,name,polse_constraints,'observability')
            ])
    eph.compute(dates)
    return eph.data.loc[ dates ]['observability'].to_list()

def get_set_rise(name,dates):
    eph = ephemeris(
        polse,
        name,
        [   
            rise_time(polse,name,'rise'),
            set_time(polse,name,'set'),
            observability(polse,name,polse_constraints,'observability')
            ])
    eph.compute(dates)
    return eph.data.loc[ dates ]['rise'].to_list(),eph.data.loc[ dates ]['set'].to_list()
    

@cli.command('rise-set')
@click.argument('name',required=True,nargs=1)
def show_rise_set(name):
    dates = pandas.date_range(start='2020-01-01',end='2020-12-31')
    
    # if not cached compute, otherwise load from cache
    sunR,sunS = get_set_rise('Sun',dates)
    R,S = get_set_rise(name,dates)

    fig=plt.figure(figsize=(20,5))
    ax=fig.add_axes([.1,.1,.8,.8])
    ax.set_ylim((0,24))
    ax.set_xlim(( dates[0], dates[-1]  ))
    ax.fill_between(dates,sunS,sunR,facecolor='black',alpha=0.5,label='night')
    ax.plot(dates,R,'o',color='blue',label='rise')
    ax.plot(dates,S,'o',color='red',label='set')
    ax.legend()
    ax.set_yticks(list(range(0,25)))
    ax.set_yticklabels(
        ['noon']+list(range(13,24))+['midnight']+list(range(1,12))+['noon'])
    months_ticks=[ datetime.datetime(2020,i,1)  for i in range(1,13)  ]
    ax.set_xticks(months_ticks)
    ax.set_xticklabels(months_labels)
    ax.set_title('%s rise-set chart' % name)
    ax.grid()
    plt.show()

@cli.command('observability')
@click.argument('name_list',required=True,nargs=-1)
def show_observability(name_list):
    '''
    Shows the observability chart for a list of objects
    '''
    dates = pandas.date_range(start='2020-01-01',end='2020-12-31')
    
    fig=plt.figure(figsize=(20,5))
    ax=fig.add_axes([.1,.1,.8,.8])
    ax.set_xlim(( dates[0],dates[-1]  ))
    
    for name in name_list:
        up_times = get_observability(name,dates)
        ax.plot(
            dates, 
            up_times,
            label=name)
    months_ticks=[ datetime.datetime(2020,i,1)  for i in range(1,13)  ]
    ax.set_xticks(months_ticks)
    ax.set_xticklabels(months_labels)
    ax.set_title("Observability chart")
    ax.grid()
    ax.legend()
    ax.set_ylabel('hours')
    plt.show()

@cli.command('sky-plot')
@click.argument('name_list',required=True,nargs=-1)
def sky_plot(name_list):
    '''
    shows the sky chart tonight at midnight
    '''
    dt = datetime.datetime(2020,12,1,0,0)
    fig=plt.figure(figsize=(5,5))
    ax=fig.add_axes([.1,.1,.8,.8],projection='polar')
    for name in name_list:
        tz = polse.timezone
        t = astropy.time.Time(tz.localize(dt))
        target = target_by_name(name,t)
        ax = astroplan.plots.plot_sky(target=target,observer=polse,time=t,ax=ax)
    plt.legend()
    plt.show()

@cli.command('wup',help='List of objects up for a given date.')
@click.option('-d','--date','date_str',required=False,help='date for which to compute wup.',
    default=str(datetime.datetime.now().date()))
def wup(date_str):
    date = datetime.datetime.strptime(date_str,'%Y-%m-%d')
    print(catalog.data)
    for name in catalog.dataindex:
        print('getting %s' % name )
        up_time = get_observability(name,[date])[0]
        if up_time > 0:
            print(name)
        

@cli.command('download',help='Download the catalogs')
def download_catalogs():
    astroquery.vizier.Vizier.ROW_LIMIT=-1
    
    # star catalog
    cat = astroquery.vizier.Vizier.get_catalogs(['I/239/hip_main'])[0].to_pandas()
    cat.to_csv('hipparcos.csv',index=False) 
    
    # galaxy catalog
    cat = astroquery.vizier.Vizier.get_catalogs(['VII/98A/catalog'])[0].to_pandas()
    cat.to_csv('VII-98A-catalog.csv',index=False) 
    
    # NGC catalog, deep-sky objects: galaxies, star clusters, and nebulae
    cat = astroquery.vizier.Vizier.get_catalogs(['VII/118/ngc2000'])[0].to_pandas()
    cat.to_csv('ngc.csv',index=False) 
    
if __name__=='__main__':
    cli()
    # example: show the observability of messier objects
    #obj_list = pandas.read_csv('object_list.csv')
    #messier_list = obj_list[ obj_list.type.eq('messier') & obj_list.magnitude.le(7)  ]
    #print(messier_list)
