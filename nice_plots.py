import matplotlib as mpl
import matplotlib.pylab as plt
from matplotlib.widgets import Slider
#import pyPLUTO as plp
#import PhyConst as phc
#import Jet_Analysis as jana
#import disk_analysis as da
import numpy as np
#import idlsave as idls # This is just for the mass_flux_plot. 
import os

from scipy.interpolate import interp1d
from scipy.ndimage import map_coordinates

### LINES
# See http://matplotlib.sourceforge.net/api/artist_api.html#module-matplotlib.lines for more
# information on line properties.
#lines.linewidth   : 2.0     # line width in points
#lines.linestyle   : -       # solid line
#lines.color       : blue
#lines.marker      : None    # the default marker
#lines.markeredgewidth  : 0.5     # the line width around the marker symbol
#lines.markersize  : 6            # markersize, in points
#lines.dash_joinstyle : miter        # miter|round|bevel
#lines.dash_capstyle : butt          # butt|round|projecting
#lines.solid_joinstyle : miter       # miter|round|bevel
#lines.solid_capstyle : projecting   # butt|round|projecting
#lines.antialiased : True         # render lines in antialised (no jaggies)



def nice_plots(fs=None):
    # Font and text properties
    if fs==None:
        mpl.rcParams['font.size'] = 24
    else:
        mpl.rcParams['font.size'] = fs
    
    mpl.rcParams['font.family']= 'serif'
    mpl.rcParams['font.serif']= 'Times'
    mpl.rcParams['text.usetex']= True
    # Figure properties
    mpl.rcParams['figure.figsize']=14,12
    mpl.rcParams['figure.subplot.left'] = .1
    mpl.rcParams['figure.subplot.right'] = .9
    mpl.rcParams['figure.subplot.bottom'] = .1
    mpl.rcParams['figure.subplot.top'] = .9
    # Lines Properties
    mpl.rcParams['lines.linewidth'] = 2.5
    mpl.rcParams['lines.markersize']  = 10
    mpl.rcParams['lines.markeredgewidth']=1.5
    # Axes and ticks properties [HINT : for getting minor ticks use plt.minorticks_on()]
    mpl.rcParams['axes.linewidth']= 3.0
    mpl.rcParams['xtick.major.size']=13.0
    mpl.rcParams['ytick.major.size']=13.0
    mpl.rcParams['xtick.minor.size']=7.0
    mpl.rcParams['ytick.minor.size']=7.0
    mpl.rcParams['xtick.major.pad']='8'
    mpl.rcParams['ytick.major.pad']='8'


def compare_images(*args,**kwargs):
    I = plp.Image()
    fig = plt.figure()
    fig.add_subplot(121)
    cdir = os.getcwd()+'/'
    if kwargs.get('pureimages',False)==True:
        var1 = args[1]
        var2 = args[2]
        Data1 = args[0]
        Data2 = args[0]
    else:    
        if kwargs.get('log',False)==True:
            var1 = np.log(args[0].__getattribute__(args[2]))
            var2 = np.log(args[1].__getattribute__(args[2]))
            Data1 = args[0]
            Data2 = args[1]
      
        else:
            var1 = args[0].__getattribute__(args[2])
            var2 = args[1].__getattribute__(args[2])
            Data1 = args[0]
            Data2 = args[1]

    if kwargs.get('scalefact')==None:
        scaling = 1.0
    else:
        scaling = kwargs.get('scalefact')
    
    I.pldisplay(scaling*var1,x1 = Data1.x1,x2=Data1.x2,cbar=(True,'vertical'))
    plt.axis([0.0,52.0,0.0,152.0])
    if kwargs.get('flines',(False,[3.0,6.0,9.0],[0.01,0.01,0.01]))[0] == True:
        I.myfieldlines(args[0],kwargs.get('flines')[1],kwargs.get('flines')[2],stream=True,colors='r',ls='-',lw=2.0)
    if kwargs.get('MHDflines',(False,319,cdir))[0] == True:
        data1MHD = plp.pload(kwargs.get('MHDflines',(True,319,cdir))[1],w_dir=kwargs.get('MHDflines',(True,319,cdir))[2])
        I.myfieldlines(data1MHD,kwargs.get('flines')[1],kwargs.get('flines')[2],stream=True,colors='w',ls='--',lw=2.0)
    plt.xlabel(kwargs.get('label1',2*['Xlabel'])[0])
    plt.ylabel(kwargs.get('label2',2*['Ylabel'])[0])
    plt.title(kwargs.get('title',2*['Title'])[0])

    fig.add_subplot(122)
    I.pldisplay(scaling*var2,x1 = Data2.x1,x2=Data2.x2,cbar=(True,'vertical'))
    plt.axis([0.0,52.0,0.0,152.0])
    if kwargs.get('flines',(False,[3.0,6.0,9.0],[0.01,0.01,0.01]))[0] == True:
        I.myfieldlines(args[1],kwargs.get('flines')[1],kwargs.get('flines')[2],stream=True,colors='r',ls='-',lw=2.0)
    if kwargs.get('MHDflines',(False,319,cdir))[0] == True:
        data2MHD = plp.pload(kwargs.get('MHDflines',(True,319,cdir))[1],w_dir=kwargs.get('MHDflines',(True,319,cdir))[2])
        I.myfieldlines(data2MHD,kwargs.get('flines')[1],kwargs.get('flines')[2],stream=True,colors='w',ls='--',lw=2.0)
    plt.xlabel(kwargs.get('label1',2*['Xlabel'])[1])
    plt.ylabel(kwargs.get('label2',2*['Ylabel'])[1])
    plt.title(kwargs.get('title',2*['Title'])[1])


def compare_plots(*args,**kwargs):
    # Arg0 - List of Data
    # Arg1 - List of Tuples [ ('variable', cutvalue, scalefact, 'log'/'semilogx'/'semilogy'/'nolog') ]
    Data_List = args[0]
    Plot_List = args[1]
    try:
        len(Data_List) == len(Plot_List)
    except:
        print 'Error Lengths should match'

    Mstar = kwargs.get('Mstar',1.0)
    ul =  kwargs.get('ul',1.0)
    urho = kwargs.get('urho',1.0)
    
    cut = kwargs.get('cut','zcut')
    Ncol = kwargs.get('Ncol', 1)

    lincol = kwargs.get('lincol','k')
    linsty = kwargs.get('linsty','-')

    pltype = kwargs.get('pltype','multiple')
    
    Vkep = np.sqrt((phc.G*phc.Msun*Mstar)/(phc.au*ul))
    Vkep = Vkep*1.0e-5 # in km/s.

    xvar=[]
    yvar=[]
   # f1 = plt.figure()
    
    

    for i in range(len(Data_List)):
        phx1 = ul*Data_List[i].x1
        phx2 = ul*Data_List[i].x2
        if cut == 'zcut':
            if Plot_List[i][0] == 'Temp':
                xvar.append(phx1[25:])
                yvar.append(Data_List[i].pr[25:,Plot_List[i][1]]/(Data_List[i].rho[25:,Plot_List[i][1]]))
            elif Plot_List[i][0] == 'Btot':
                xvar.append(phx1)
                Bt = np.sqrt(Data_List[i].b1**2 +Data_List[i].b2**2 + Data_List[i].b3**2)  
                yvar.append(Bt[:,Plot_List[i][1]])
            else:
                xvar.append(phx1)
                yvar.append(Data_List[i].__getattribute__(Plot_List[i][0])[:,Plot_List[i][1]])
        else:
            xvar.append(phx2)
            if Plot_List[i][0] == 'Temp':
                yvar.append(Data_List[i].pr[Plot_List[i][1],:]/(Data_List[i].rho[Plot_List[i][1],:]))
            elif Plot_List[i][0] == 'Btot':
               Bt = np.sqrt(Data_List[i].b1**2 +Data_List[i].b2**2 + Data_List[i].b3**2)  
               yvar.append(Bt[Plot_List[i][1],:])
            else:
                yvar.append(Data_List[i].__getattribute__(Plot_List[i][0])[Plot_List[i][1],:])

    for i in range(len(xvar)):
        
        if pltype == 'single':
            plt.subplot(111)
        else:
            plt.subplot(len(xvar),Ncol,i+1)
        

        if Plot_List[i][3]=='nolog':
            plt.plot(xvar[i],Plot_List[i][2]*yvar[i],color=lincol[i],linestyle=linsty[i])
        elif Plot_List[i][3]=='semilogx':
            semilogx(xvar[i],Plot_List[i][2]*yvar[i],color=lincol[i],linestyle=linsty[i])
        elif Plot_List[i][3]=='semilogy':
            plt.semilogy(xvar[i],Plot_List[i][2]*yvar[i],color=lincol[i],linestyle=linsty[i])
        else:
            plt.loglog(xvar[i],Plot_List[i][2]*yvar[i],color=lincol[i],linestyle=linsty[i])

            
            
        
        
        


def symimage(*args,**kwargs):
    I = plp.Image()
    T = plp.Tools()
    #fig = plt.figure()
    Dmhd = plp.pload(319,w_dir='/Users/bhargavvaidya/pyPLUTO-wdir/M30a055b5_data/')
    Data = args[0]
    sclf =  kwargs.get('scalefact',1.0)
    if kwargs.get('log',False)==True:
        var = np.log10(sclf*Data.__getattribute__(args[1]))
    else:
        var = sclf*Data.__getattribute__(args[1])
    x1 = Data.x1
    x2 = Data.x2                   
    plt.axis([-x1.max(),x1.max(),0.0,x2.max()])

    plt.title(kwargs.get('title',"Title"),size=kwargs.get('size'))
    plt.xlabel(kwargs.get('label1',"Xlabel"),size=kwargs.get('size'),labelpad=6)
    plt.ylabel(kwargs.get('label2',"Ylabel"),size=kwargs.get('size'),labelpad=6)
    
    xm,ym = np.meshgrid(x1.T,x2.T)
    xmc = T.congrid(xm,2*([kwargs.get('arrows')[1]]),method='linear')
    ymc = T.congrid(ym,2*([kwargs.get('arrows')[1]]),method='linear')
    v1c = T.congrid(Data.v1.T,2*([kwargs.get('arrows')[1]]),method='linear')
    v2c = T.congrid(Data.v2.T,2*([kwargs.get('arrows')[1]]),method='linear')
    v1cn,v2cn = v1c/np.sqrt(v1c**2+v2c**2),v2c/np.sqrt(v1c**2+v2c**2)
                        
    
    if kwargs.get('ysym',False) == True:
        ysymvar = np.flipud(var)
        plc1=plt.pcolormesh(-x1[::-1],x2,ysymvar.T,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        plc2=plt.pcolormesh(x1,x2,var.T,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        
        if kwargs.get('cbar',(False,''))[0] == True:
            plcol=plt.colorbar(orientation=kwargs.get('cbar')[1])
        if kwargs.get('arrows',(False,20))[0]==True:
            plq1=plt.quiver(xmc,ymc,v1cn,v2cn,color='k')
            plq2=plt.quiver(-xmc,ymc,-v1cn,v2cn,color='k')
            

        if kwargs.get('flines',(False,[3.0,6.0,10.0,15.0],[0.001,0.001,0.001,0.001]))[0]==True:
            x0arr = kwargs.get('flines')[1]
            y0arr = kwargs.get('flines')[2]
            for i in range(np.size(x0arr)):
                flmhd = I.field_line(Dmhd.b1,Dmhd.b2,Dmhd.x1,Dmhd.x2,Dmhd.dx1,Dmhd.dx2,x0arr[i],y0arr[i])
                Qxmhd = flmhd.get('qx')
                Qymhd = flmhd.get('qy')
                fl =I.field_line(Data.b1,Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx2,x0arr[i],y0arr[i])
                Qx= fl.get('qx')
                Qy= fl.get('qy')
                plfl1=plt.plot(Qx,Qy,color=kwargs.get('colors','r'),ls=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
                plfl2=plt.plot(-Qx,Qy,color=kwargs.get('colors','r'),ls=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
                plfl3=plt.plot(Qxmhd,Qymhd,color=kwargs.get('colors','w'),ls=kwargs.get('ls','--'),lw=kwargs.get('lw',2.0))
                plfl4=plt.plot(-Qxmhd,Qymhd,color=kwargs.get('colors','w'),ls=kwargs.get('ls','--'),lw=kwargs.get('lw',2.0))

    if kwargs.get('xsym',False) == True:
        xsymvar = np.fliplr(var)
        plc1=plt.pcolormesh(x1,-x2[::-1],xsymvar.T,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        plc2=plt.pcolormesh(x1,x2,var.T,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        if kwargs.get('cbar',(False,''))[0] == True:
            plcol=plt.colorbar(orientation=kwargs.get('cbar')[1])

        if kwargs.get('arrows',(False,20))[0]==True:
            plq1=plt.quiver(xmc,ymc,v1cn,v2cn,color='g')
            plq2=plt.quiver(xmc,-ymc,v1cn,-v2cn,color='g')
            
            
        if kwargs.get('flines',(False,[3.0,6.0,10.0,15.0],[0.001,0.001,0.001,0.001]))[0]==True:
            x0arr = kwargs.get('flines')[1]
            y0arr = kwargs.get('flines')[2]
            for i in range(np.size(x0arr)):
                flmhd = I.field_line(Dmhd.b1,Dmhd.b2,Dmhd.x1,Dmhd.x2,Dmhd.dx1,Dmhd.dx2,x0arr[i],y0arr[i])
                Qxmhd = flmhd.get('qx')
                Qymhd = flmhd.get('qy') 
                fl =I.field_line(Data.b1,Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx2,x0arr[i],y0arr[i])
                Qx= fl.get('qx')
                Qy= fl.get('qy')
                plfl1= plt.plot(Qx,Qy,color=kwargs.get('colors','r'),ls=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
                plfl2= plt.plot(Qx,-Qy,color=kwargs.get('colors','r'),ls=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
                plfl1= plt.plot(Qxmhd,Qymhd,color=kwargs.get('colors','w'),ls=kwargs.get('ls','--'),lw=kwargs.get('lw',2.0))
                plfl2= plt.plot(Qxmhd,-Qymhd,color=kwargs.get('colors','w'),ls=kwargs.get('ls','--'),lw=kwargs.get('lw',2.0))
                
                

    if kwargs.get('fullsym',False) == True:
        ysymvar = np.flipud(var)
        plc1=plt.pcolormesh(-x1[::-1],x2,ysymvar.T,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        xsymvar = np.fliplr(var)
        plc2=plt.pcolormesh(x1,-x2[::-1],xsymvar.T,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        plc3=plt.pcolormesh(-x1,-x2,var.T,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        plc4=plt.pcolormesh(x1,x2,var.T,vmin=kwargs.get('vmin',np.min(var)),vmax=kwargs.get('vmax',np.max(var)))
        if kwargs.get('cbar',(False,''))[0] == True:
            plcol=plt.colorbar(orientation=kwargs.get('cbar')[1])
            
        if kwargs.get('arrows',(False,20))[0]==True:
            plq1=plt.quiver(-xmc,ymc,-v1cn,v2cn,color='w')
            plq2=plt.quiver(xmc,-ymc,v1cn,-v2cn,color='w')
            plq3=plt.quiver(xmc,ymc,v1cn,v2cn,color='w')
            plq4=plt.quiver(-xmc,-ymc,-v1cn,-v2cn,color='w')
            
        if kwargs.get('flines',(False,[2.0,4.0,7.0,13.0],[0.001,0.001,0.001,0.001]))[0]==True:
            x0arr = kwargs.get('flines')[1]
            y0arr = kwargs.get('flines')[2]
            for i in range(np.size(x0arr)):
                fl =I.field_line(Data.b1,Data.b2,Data.x1,Data.x2,Data.dx1,Data.dx2,x0arr[i],y0arr[i])
                Qx= fl.get('qx')
                Qy= fl.get('qy')
                plfl1=plt.plot(Qx,Qy,color=kwargs.get('colors','r'),ls=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
                plfl2=plt.plot(-Qx,Qy,color=kwargs.get('colors','r'),ls=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
                plfl3=plt.plot(Qx,-Qy,color=kwargs.get('colors','r'),ls=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
                plfl4=plt.plot(-Qx,-Qy,color=kwargs.get('colors','r'),ls=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
        
    
    return [plc1,plc2,plcol,plq1,plq2,plfl1,plfl2,plfl3,plfl4] 



def zetaper(xfl=None,nmhd=None,nrad=None,w_dir=None):
        if xfl==None :
            print "The time step for MHD run not given - Default xfl=5.0"
            xfl=5.0
        if nmhd == None :
            print "The time step for MHD run not given - Default 319"
            nmhd = 319
        if nrad == None :
            print "The time step for RAD run not given - Default 319"
            nrad = 319
	if w_dir==None :
	    print w_dir

        dMHD = plp.pload(nmhd,w_dir=w_dir)
        dRAD = plp.pload(nrad,w_dir=w_dir)
        I = plp.Image()
        
        fldictMHD = I.field_line(dMHD.b1,dMHD.b2,dMHD.x1,dMHD.x2,dMHD.dx1,dMHD.dx2,xfl,0.01)
        ZetaMHD = (180.0/np.pi)*np.arctan2(fldictMHD.get('qx'),fldictMHD.get('qy'))
        fldictRAD = I.field_line(dRAD.b1,dRAD.b2,dRAD.x1,dRAD.x2,dRAD.dx1,dRAD.dx2,xfl,0.01)
        ZetaRAD = (180.0/np.pi)*np.arctan2(fldictRAD.get('qx'),fldictRAD.get('qy'))
        nlenmin = np.min([len(ZetaRAD),len(ZetaMHD)])
        Diffinper = 100.0*((ZetaRAD[:nlenmin]-ZetaMHD[:nlenmin])/(ZetaMHD[:nlenmin]))
        return [dMHD.x2[:nlenmin].T,Diffinper.T]

def force_plots(Data,**kwargs):
    Fcl = jana.Force()
   
    var_fl=Fcl.Fline_values(Data,**kwargs)
    Res_Grav= np.sqrt(var_fl['Fl_Gr'][0,:]**2 + var_fl['Fl_Gr'][1,:]**2)
    Res_Press= np.sqrt(var_fl['Fl_Pr'][0,:]**2 + var_fl['Fl_Pr'][1,:]**2)
    Res_Centr = np.sqrt(var_fl['Fl_Cf'][0,:]**2 + var_fl['Fl_Cf'][1,:]**2)
    Res_Lor = np.sqrt(var_fl['Fl_Lf'][0,:]**2 + var_fl['Fl_Lf'][1,:]**2)
    Res_StRad = np.sqrt(var_fl['Fl_StRf'][0,:]**2 + var_fl['Fl_StRf'][1,:]**2)
    tfl_Grav = var_fl['tFl_Gr'][0,:]
    tfl_Press = var_fl['tFl_Pr'][0,:]
    tfl_Centri = var_fl['tFl_Cf'][0,:]
    tfl_bpolpr = var_fl['tFl_Lf1'][0,:]
    tfl_bphipr = var_fl['tFl_Lf2'][0,:]
    tfl_pinch = var_fl['tFl_Pinch'][0,:]
    tfl_StRf = var_fl['tFl_StRf'][0,:]

    afl_Grav = var_fl['aFl_Gr'][0,:]
    afl_Press = var_fl['aFl_Pr'][0,:]
    afl_Centri = var_fl['aFl_Cf'][0,:]
    afl_bpolpr = var_fl['aFl_Lf1'][0,:]
    afl_bphipr = var_fl['aFl_Lf2'][0,:]
    afl_pinch = var_fl['aFl_Pinch'][0,:]
    afl_StRf = var_fl['aFl_StRf'][0,:]
    
    
    if kwargs.get('DiskRad',False)==True:
        Res_DiskRad = np.sqrt(var_fl['Fl_DiskRf'][0,:]**2 + var_fl['Fl_DiskRf'][1,:]**2)
   
    Qy = var_fl['Qy']

    if kwargs.get('Gravity',False)==True:
        plt.loglog(Qy,Res_Grav,color=kwargs.get('color','black'),linestyle=kwargs.get('ls','solid'))

    if kwargs.get('Centrifugal',False)==True:
        plt.loglog(Qy,Res_Centr,color=kwargs.get('color','black'),linestyle=kwargs.get('ls','solid'))

    if kwargs.get('Pressure',False)==True:
        plt.loglog(Qy,Res_Press,color=kwargs.get('color','black'),linestyle=kwargs.get('ls','solid'))

    if kwargs.get('Lorentz',False)==True:
        plt.loglog(Qy,Res_Lor,color=kwargs.get('color','black'),linestyle=kwargs.get('ls','solid'))

    if kwargs.get('StRad',False)==True:
        plt.loglog(Qy,Res_StRad,color=kwargs.get('color','black'),linestyle=kwargs.get('ls','solid'))

    if kwargs.get('DiskRad',False)==True:
        plt.loglog(Qy,Res_DiskRad,color=kwargs.get('color','black'),linestyle=kwargs.get('ls','solid'))

    if kwargs.get('Transf',False)==True:
        pgr = plt.loglog(Qy,tfl_Grav+np.abs(tfl_pinch),color='green',linestyle=kwargs.get('ls','solid'))
       # pstr = plt.loglog(Qy,np.abs(tfl_StRf),color='blue',linestyle=kwargs.get('ls','solid'))
       # ppr = plt.loglog(Qy,np.abs(tfl_Press),color='black',linestyle=kwargs.get('ls','solid'))
        pcf = plt.loglog(Qy,tfl_Centri+np.abs(tfl_Press)+np.abs(tfl_StRf)+np.abs(tfl_bphipr),color='red',linestyle=kwargs.get('ls','solid'))
        #pbpolpr = plt.loglog(Qy,np.abs(tfl_bpolpr),color='green',linestyle=kwargs.get('ls','solid'))
        #pbphipr = plt.loglog(Qy,np.abs(tfl_bphipr),color='blue',linestyle=kwargs.get('ls','solid'))
        #ppinch = plt.loglog(Qy,np.abs(tfl_pinch),color='magenta',linestyle=kwargs.get('ls','solid'))

        #plt.legend([pbpolpr,pbphipr,ppinch],[r'$\nabla_{\perp}(B_{p}^{2}/2)$',r'$\nabla_{\perp}(B_{\phi}^{2}/2)$',r'$(B_{\phi}^{2}/r) \nabla_{\perp} r$'],'lower left')
        #plt.legend([pgr,ppr,pcf,pstr],[r'$\rho \nabla_{\perp}\phi$',r'$\nabla_{\perp}P$',r'$\rho\Omega^{2}r \nabla_{\perp}r$',r'$\nabla_{\perp}(B^{2}/2)$',r'$(B_{\phi}^{2}/r) \nabla_{\perp} r$'],'lower left')

    if kwargs.get('Paraf',False)==True:
        pgr = plt.loglog(Qy,afl_Grav,color='green',linestyle=kwargs.get('ls','solid'))
        #pstr = plt.loglog(Qy,np.abs(afl_StRf),color='brown',linestyle=kwargs.get('ls','solid'))
        ppr = plt.loglog(Qy,np.abs(afl_Press),color='black',linestyle=kwargs.get('ls','solid'))
        pcf = plt.loglog(Qy,afl_Centri,color='red',linestyle=kwargs.get('ls','solid'))
        pbpolpr = plt.loglog(Qy,np.abs(afl_bpolpr),color='cyan',linestyle=kwargs.get('ls','solid'))
        pbphipr = plt.loglog(Qy,np.abs(afl_bphipr),color='blue',linestyle=kwargs.get('ls','solid'))
        ppinch = plt.loglog(Qy,np.abs(afl_pinch),color='magenta',linestyle=kwargs.get('ls','solid'))

        
        

    

def ldinstable(nparr=None,nmhd=None,w_dir=None,**kwargs):
    if nmhd == None : nmhd = 319
    if nparr == None : nparr = [198,199,200]
    Vref = 1.0e-5*np.sqrt((phc.G*phc.Msun*kwargs.get('Mstar',30.0))/(phc.au*kwargs.get('ul',1.0)))
    print '---------------------------------------------'
    print 'Mass of Star :',kwargs.get('Mstar',30)
    print 'Inner disk radius :',kwargs.get('ul',1.0)
    print 'UNIT VELOCITY [km/s] :',Vref
    print '----------------------------------------------'
  


    Fc = jana.Force()
    Data_List=[plp.pload(nmhd,w_dir=w_dir)]
    for j in nparr:
        Data_List.append(plp.pload(j,w_dir=w_dir))

    Vp_Dict_List =[]
    for p in Data_List:
        Vp_Dict_List.append(Fc.proj_force(p,p.v1,p.v2,**kwargs))

    f1 = plt.figure(num=1)
    ax2 = f1.add_subplot(212)
    plt.plot( Vp_Dict_List[0]['Qy'], Vref*Vp_Dict_List[0]['para_flvalues'][0,:],'k--')
    for q in  Vp_Dict_List[1:]:
        plt.plot(q['Qy'],Vref*q['para_flvalues'][0,:])
    
    plt.xlabel(r'Distance Along the Field Line : s [AU]',labelpad=6)
    plt.ylabel(r'Poloidal Velocity : $V_{\rm pol} [\rm km\, \rm s^{-1}]$',labelpad=10)

    
    
        
                
        

def deg_coll_plots(dirlist=None, **kwargs):
    if kwargs.get("perzeta",False)==True:
        zetalist = []
        for dirnames in dirlist:
          zetalist.append(zetaper(w_dir=dirnames,nrad=kwargs.get('nrad',638)))
          #f1 = plt.figure()
        plt.xlabel(r'Distance Along Field Line : s [AU]',labelpad=6)
        plt.ylabel(r'$\Delta \zeta [\%]$',labelpad=10)
        plt.plot(zetalist[0][0],zetalist[0][1],'k-.')
        plt.plot(zetalist[1][0],zetalist[1][1],'k--')
        plt.plot(zetalist[2][0],zetalist[2][1],'k')
        plt.axis([0.0,150.0,0.0,40.0])
        plt.minorticks_on()
    else:
        magalf = np.linspace(19.67,19.67,100)
        magfast = np.linspace(12.53,12.53,100)
        
        Msarr = [20.0,25.0,30.0,45.0,50.0,60.0]
        alfMs = [20.44,23.27,26.05, 26.84, 28.86,31.45]
        fastMs = [15.57,19.94,23.40,24.35,26.03,29.25]
        asymMs=[]
        Alparr = [0.55,0.60,0.65]
        alfAlp =[31.35,23.73,20.93]
        fastAlp=[29.25,21.48,17.28]
        asymAlp =[]
        Denarr =[3.0e-14,5.0e-14,1.0e-13,5.0e-13]
        alfDen =[30.30,26.05,22.91,20.42]
        fastDen =[27.79,23.40,19.47,15.03]
        asymDen=[]
        Betarr=[1.0,3.0,5.0]
        alfBet =[21.98,25.63,26.05]
        fastBet=[16.45,20.68,23.40]
        asymBet=[]

        #nice_plots()
        #f1 = plt.figure()
        #ax1 = f1.add_subplot(111)
        
        plt.ylabel(r"Opening Angle : $\phi^\circ$",labelpad=10)
        plt.minorticks_on()
        if kwargs.get('cplot',False)==True:
            dum1 = np.linspace(1.0e-14,1.0e-12,100)
            ms1 = plt.semilogx(Denarr,alfDen,'r*')
            ms2 = plt.semilogx(Denarr,fastDen,'bo')
            plt.xlabel(r"Base Density : $\rho_0$ [g/cm$^{3}$]")
            plt.semilogx(dum1,magalf,'r--',dum1,magfast,'b--')
            plt.axis([1.0e-14,1.0e-12,10.0,35.0])
           # plt.legend([ms1,ms2],[r'Alfven point',r'Fast point'],loc='lower left')
        else:
            dum1 = np.linspace(10.0,70.0,100)
            ms1 = plt.plot(Msarr,alfMs,'k*')
            ms2 = plt.plot(Msarr,fastMs,'ko')
            plt.xlabel(r"Stellar Mass : $M_{*}$ [$M_{\odot}$]",labelpad=6)
            plt.plot(dum1,magalf,'k--',dum1,magfast,'k')
            plt.axis([10.0,70.0,10.0,35.0])

            
            ## f2 = plt.figure()
            ## ax1 = f2.add_subplot(211)
            #dumalf = np.linspace(0.5,0.7,100)
            #ms1 = plt.plot(Alparr,alfAlp,'k*')
            #ms2 = plt.plot(Alparr,fastAlp,'ko')
            #plt.xlabel(r"Line Force Parameter : $\alpha$",labelpad=6)
            #plt.ylabel(r"Opening Angle : $\phi^\circ$",labelpad=10)
            #plt.minorticks_on()
            #plt.plot(dumalf,magalf,'k--',dumalf,magfast,'k')
            #plt.axis([0.5,0.7,10.0,35.0])
            
            ## ax2 = f2.add_subplot(212)
            #dumrho = np.linspace(1.0e-14,1.0e-12,100)
            #ms1 = plt.semilogx(Denarr,alfDen,'k*')
            #ms2 = plt.semilogx(Denarr,fastDen,'ko')
            #plt.xlabel(r"Base Density : $\rho_0$ [g cm$^{-3}$]",labelpad=6)
            #plt.ylabel(r"Opening Angle : $\phi^\circ$",labelpad=10)
            #plt.minorticks_on()
            #plt.semilogx(dumrho,magalf,'k--',dumrho,magfast,'k')
            #plt.axis([1.0e-14,1.0e-12,10.0,35.0])

            
            
            


def conto_plots(Data,DataMHD,**kwargs):
    Qu = jana.quantities()
    Cur = Qu.Current(Data)
    CurMHD = Qu.Current(DataMHD)
    Fastsp = Qu.Magspeed(Data)['fast']
    Slowsp = Qu.Magspeed(Data)['slow']
    Alfvsp = Qu.Magspeed(Data)['alfven']
    Vpol = np.sqrt(Data.v1**2 + Data.v2**2)

    FastspMHD = Qu.Magspeed(DataMHD)['fast']
    SlowspMHD = Qu.Magspeed(DataMHD)['slow']
    AlfvspMHD = Qu.Magspeed(DataMHD)['alfven']
    VpolMHD = np.sqrt(DataMHD.v1**2 + DataMHD.v2**2)

    fastrat = Fastsp/Vpol
    slowrat = Slowsp/Vpol
    Alfvrat = Alfvsp/Vpol

    fastratMHD = FastspMHD/VpolMHD
    slowratMHD = SlowspMHD/VpolMHD
    AlfvratMHD = AlfvspMHD/VpolMHD

    

    f1 = plt.figure(num=1)
    ax1 = f1.add_subplot(121)
    currcont=plt.contour(Data.x1,Data.x2,Cur.T,kwargs.get('Currents',[-0.7,-0.6,-0.5,-0.4,-0.3]),colors=kwargs.get('colors','r'),linestyles=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
    plt.clabel(currcont,manual=True)
    currmhdcont = plt.contour(DataMHD.x1,DataMHD.x2,CurMHD.T,kwargs.get('Currents',[-0.7,-0.6,-0.5,-0.4,-0.3]),colors=kwargs.get('colors','k'),linestyles=kwargs.get('ls','-'),lw=kwargs.get('lw',2.0))
    plt.clabel(currmhdcont,manual=True)
    plt.xlabel(r'r [AU]')
    plt.ylabel(r'z [AU]')
    plt.title(r'Total Current -- $\int\int J_{z} r dr d\phi$')
    ax2 = f1.add_subplot(122)
    plt.xlabel(r'r [AU]')
    #plt.ylabel(r'$z [AU]$')
    plt.title(r'Critical Surfaces')
    
    fcont=plt.contour(Data.x1,Data.x2,fastrat.T,[1.0],colors='r',linestyles='solid')
    plt.clabel(fcont,inline=1,fmt=r'Fast')
    scont=plt.contour(Data.x1,Data.x2,slowrat.T,[1.0],colors='r',linestyles='dashdot')
    plt.clabel(scont,inline=1,fmt=r'Slow')
    acont=plt.contour(Data.x1,Data.x2,Alfvrat.T,[1.0],colors='r',linestyles='dashed')
    plt.clabel(acont,inline=1,fmt=r'Alfv$\acute{e}$n')

    mfcont=plt.contour(DataMHD.x1,DataMHD.x2,fastratMHD.T,[1.0],colors='k',linestyles='solid')
    plt.clabel(mfcont,manual=1,fmt=r'Fast')
    mscont=plt.contour(DataMHD.x1,DataMHD.x2,slowratMHD.T,[1.0],colors='k',linestyles='dashdot')
    plt.clabel(mscont,manual=1,fmt=r'Slow')
    macont=plt.contour(DataMHD.x1,DataMHD.x2,AlfvratMHD.T,[1.0],colors='k',linestyles='dashed')
    plt.clabel(macont,manual=1,fmt=r'Alfv$\acute{e}$n')

    
def mass_flux_plot(*args,**kwargs):
    fltm = idls.read(args[0])
    injm = idls.read(args[1])
    f1 = plt.figure()

    ax1 = f1.add_subplot(211)
    plt.plot(injm.nt_sc,injm.nmf_rscale,'r')
    plt.plot(injm.nt_sc,injm.nmf_zscale,'b')
    plt.plot(injm.nt_sc,injm.nmf_z0scale,'k')
    plt.plot(injm.nt_sc,(injm.nmf_rscale+injm.nmf_zscale),'g')
    plt.axis([0.0,160.0,0.0,3.5e-5])
    plt.minorticks_on()
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.1f" % x, locs*1e5))
    plt.text(0.0, 1.03, r'$10^{-5}$', transform = plt.gca().transAxes)
    plt.xlabel(r'Time [yr]',labelpad=6)
    plt.ylabel(r'$\dot{\rm M}_{\rm out} [ \rm{M}_{\odot} \rm{yr}^{-1} ]$',labelpad=15)
    
    ax2 = f1.add_subplot(212)
    plt.plot(fltm.nt_sc,fltm.nmf_rscale,'r')
    plt.plot(fltm.nt_sc,fltm.nmf_zscale,'b')
    plt.plot(fltm.nt_sc,fltm.nmf_z0scale,'k')
    plt.plot(fltm.nt_sc,(fltm.nmf_rscale+fltm.nmf_zscale),'g')
    plt.axis([0.0,160.0,0.0,4.0e-5])
    plt.minorticks_on()
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.1f" % x, locs*1e5))
    plt.text(0.0, 1.03, r'$10^{-5}$', transform = plt.gca().transAxes)
    plt.xlabel(r'Time [yr]',labelpad=6)
    plt.ylabel(r'$\dot{\rm M}_{\rm out} [ \rm {M}_{\odot} \rm{yr}^{-1} ]$',labelpad=15)

def minj_mflux(*args,**kwargs):
    Qu = jana.quantities()
    Minj_List =[]
    Minj_MHD_List=[]
    for i in range(len(args[1])):
        Minj_List.append(Qu.Mflux(args[1][i],Mstar=args[0][i],scale=True)['Mfr']+Qu.Mflux(args[1][i],Mstar=args[0][i],scale=True)['Mfz'])
        
    MHD_30minj = Qu.Mflux(args[2],Mstar=30.0,scale=True)['Mfr']+Qu.Mflux(args[2],Mstar=30.0,scale=True)['Mfz']
    for i in args[0]:
        Minj_MHD_List.append(MHD_30minj*np.sqrt(i/30.0))
    
    f1 = plt.figure(num=1)
    ax1 = f1.add_subplot(211)
    plt.axis([10.0,70.0,2.0e-5,7.99e-5])
    plt.plot(args[0],Minj_MHD_List,'k*')
    plt.plot(args[0],Minj_List,'ko')
    plt.minorticks_on()
    locs,labels = plt.yticks()
    plt.yticks(locs, map(lambda x: "%.1f" % x, locs*1e5))
    plt.text(0.0, 1.03, r'$10^{-5}$', transform = plt.gca().transAxes)
    plt.xlabel(r'Stellar Mass : $M_{*} [M_{\odot}]$')
    plt.ylabel(r' $\dot{M}_{\rm vert} + \dot{M}_{\rm rad}\, [M_{\odot}\,\rm yr^{-1}]$')
    
    ax2 = f1.add_subplot(212)
    plt.axis([10.0,70.0,0.0,50.0])
    plt.plot(args[0],100*((np.array(Minj_List)-np.array(Minj_MHD_List))/np.array(Minj_MHD_List)),'k^')
    plt.minorticks_on()
    plt.xlabel(r'Stellar Mass : $M_{*} [M_{\odot}]$')
    plt.ylabel(r'$\%$ Change in Total Mass Outflow Rates')
    
    

def para_per_plot(*args,**kwargs):
    #Getting the parallel component of the force.
    Fc = jana.Force()
    D = args[0]
    Gdict=Fc.Gravity(D)
    Pdict=Fc.Pressure(D)
    Cdict=Fc.Centrifugal(D)
    Ldict=Fc.Lorentz(D)
    MagPdict = Fc.Mag_Pressure(D)
    StRdict=Fc.Stellar_Rad(D,**kwargs)
    fscale = (phc.G*phc.Msun*kwargs.get('Mstar',30)/(kwargs.get('ul',1.0)*phc.au)**2)
    print fscale

    #1. Gravity Projection
    grav_proj = Fc.proj_force(D,Gdict['G_r'],Gdict['G_z'])
    print '>> Finished Projecting Gravity force'
    #2. Pressure Projection
    press_proj = Fc.proj_force(D,Pdict['Fp_r'],Pdict['Fp_z'])
    print '>> Finished Projecting Pressure force'
    #3. Centrifugal Projection
    cf_proj = Fc.proj_force(D,Cdict['Fcf_r'],Cdict['Fcf_z'])
    print '>> Finished Projecting Centrifugal force'
    #4. Lorentz Projection
    lf_proj = Fc.proj_force(D,Ldict['Fl_r'],Ldict['Fl_z'])
    print '>> Finished Projecting Lorentz force'

    dum1 = np.linspace(1e-6,1e-0,100) 
    #dum1 = np.linspace(-0.02,1.0,100)
    fastrad = np.linspace(9.48,9.48,100)
    fastmhd = np.linspace(38.59,38.59,100)
    alfmhd = np.linspace(4.41,4.41,100)
    alfrad= np.linspace(1.71,1.71,100)

    f1 = plt.figure(num=1,figsize=(14,10))
    ax1 = f1.add_subplot(121)
    
    ax1.axis([0.5,150.0,1.0e-4,1.0e-0])
    #ax1.axis([0.5,150.0,-0.01,0.01])
    #ax1.set_aspect('equal')
    #if kwargs.get('Parallel',False)==True:
    print "GOING TO PLOT THE PARALLEL COMPONENTS"
    #plt.loglog(grav_proj['Qy'],grav_proj['para_flvalues'][0,:],'k',linestyle=kwargs.get('ls','-'))
    plt.loglog(press_proj['Qy'],fscale*press_proj['para_flvalues'][0,:],'r',linestyle=kwargs.get('ls','-'))
    plt.loglog(cf_proj['Qy'],fscale*cf_proj['para_flvalues'][0,:],'b',linestyle=kwargs.get('ls','-'))
    plt.loglog(lf_proj['Qy'],fscale*lf_proj['para_flvalues'][0,:],'g',linestyle=kwargs.get('ls','-'))
   
    #plt.plot(press_proj['Qy'],press_proj['para_flvalues'][0,:]-grav_proj['para_flvalues'][0,:],'r',linestyle=kwargs.get('ls','-'))
    #plt.plot(cf_proj['Qy'],cf_proj['para_flvalues'][0,:],'b',linestyle=kwargs.get('ls','-'))
    #plt.plot(lf_proj['Qy'],lf_proj['para_flvalues'][0,:],'g',linestyle=kwargs.get('ls','-'))
    #plt.plot(press_proj['Qy'],press_proj['para_flvalues'][0,:]-grav_proj['para_flvalues'][0,:]+cf_proj['para_flvalues'][0,:]+lf_proj['para_flvalues'][0,:],'k',linestyle=kwargs.get('ls','-'))

    if kwargs.get('StForce',False)==True:
        #5. StRad Projection
        str_proj = Fc.proj_force(D,StRdict['Fr_r'],StRdict['Fr_z'])
        print '>> Finished Projecting StRad force'
        plt.loglog(str_proj['Qy'],fscale*str_proj['para_flvalues'][0,:],'brown',linestyle=kwargs.get('ls','-'))

    plt.xlabel(r'Distance Along the Field Line : s [AU]',labelpad=6)
    plt.ylabel(r'Specific Forces Parallel to Field line [cm s$^{-2}$]',labelpad=8)
    plt.loglog(fastmhd,dum1,'m--')
    plt.loglog(alfmhd,dum1,'k--')
    plt.loglog(fastrad,dum1,'m-')
    plt.loglog(alfrad,dum1,'k-')
    #plt.plot(fastmhd,dum1,'m--')
    #plt.plot(alfmhd,dum1,'k--')
    
    
    

    ax2 = f1.add_subplot(122)
    #ax2.set_aspect('equal')
    ax2.axis([0.5,150.0,1.0e-4,1.0e-0])
    #ax2.axis([0.5,150.0,-0.02,0.02])
    #else:
    print "GOING TO PLOT THE PERPENDICULAR COMPONENTS"
    #plt.loglog(grav_proj['Qy'],grav_proj['perp_flvalues'][0,:],'k',linestyle=kwargs.get('ls','-'))
    plt.loglog(press_proj['Qy'],fscale*press_proj['perp_flvalues'][0,:],'r',linestyle=kwargs.get('ls','-'))
    plt.loglog(cf_proj['Qy'],fscale*cf_proj['perp_flvalues'][0,:],'b',linestyle=kwargs.get('ls','-'))
    plt.loglog(lf_proj['Qy'],fscale*lf_proj['perp_flvalues'][0,:],'g',linestyle=kwargs.get('ls','-'))

   # plt.plot(press_proj['Qy'],press_proj['perp_flvalues'][0,:]-grav_proj['perp_flvalues'][0,:],'r',linestyle=kwargs.get('ls','-'))
  #  plt.plot(cf_proj['Qy'],cf_proj['perp_flvalues'][0,:],'b',linestyle=kwargs.get('ls','-'))
  #  plt.plot(lf_proj['Qy'],lf_proj['perp_flvalues'][0,:],'g',linestyle=kwargs.get('ls','-'))
   # plt.plot(press_proj['Qy'],press_proj['perp_flvalues'][0,:]-grav_proj['perp_flvalues'][0,:]+cf_proj['perp_flvalues'][0,:]+lf_proj['perp_flvalues'][0,:],'k',linestyle=kwargs.get('ls','-'))


    if kwargs.get('StForce',False)==True:
        #5. StRad Projection
        str_proj = Fc.proj_force(D,StRdict['Fr_r'],StRdict['Fr_z'])
        print '>> Finished Projecting StRad force'
        plt.loglog(str_proj['Qy'],fscale*str_proj['perp_flvalues'][0,:],'brown',linestyle=kwargs.get('ls','-'))

    plt.xlabel(r'Distance Along the Field Line : s [AU]',labelpad=6)
    plt.ylabel(r'Specific Forces Perpendicular to Field line [cm s$^{-2}$]',labelpad=6)
    plt.loglog(fastmhd,dum1,'m--')
    plt.loglog(alfmhd,dum1,'k--')
    plt.loglog(fastrad,dum1,'m-')
    plt.loglog(alfrad,dum1,'k-')
    #plt.plot(fastmhd,dum1,'m--')
    #plt.plot(alfmhd,dum1,'k--')




    
    
    
    
    

    
    

def mag_forces(*args,**kwargs):
    Tool = plp.Tools()
    Fc = jana.Force()
    I = plp.Image()
    D = args[0]

    #Lorentz_Force_Dict = Fc.Lorentz(D)
    StForce_Dict = Fc.Pressure(D)

    
    Flr = StForce_Dict['Fp_r']
    Flz = StForce_Dict['Fp_z']
    
    #Flr = Lorentz_Force_Dict['Fl_r']
    #Flz = Lorentz_Force_Dict['Fl_z']
    
    spx= kwargs.get('spacingx',40)
    spy= kwargs.get('spacingx',80)
    nxr = D.n1/spx
    nyr = D.n2/spy
    nshape = ([nxr,nyr])

    b1new = Tool.congrid(D.b1,nshape) 
    b2new = Tool.congrid(D.b2,nshape)
    x1new = Tool.congrid(D.x1,([nxr]))
    x2new = Tool.congrid(D.x2,([nyr]))

    
 
    Flr_p = Tool.congrid(Flr,nshape)
    Flz_p = Tool.congrid(Flz,nshape)
    magpol_p = np.sqrt(Flr_p*Flr_p + Flz_p*Flz_p)


    #plt.quiver(x1new,x2new,(Flr_p/magpol_p).T,(Flz_p/magpol_p).T,color=kwargs.get('color','k'),headwidth=kwargs.get('hw',3))

    #FOR THE PARALLEL FORCES WE HAVE TO PROJECT THE Fl TO THE Bp
    if kwargs.get('Para',False)==True:
        phi = np.arctan2(b2new,b1new)
        f2 = plt.figure(num=2)
        
        theta = np.zeros(phi.shape)
        alpha = np.arctan2(Flz_p,Flr_p)

        

        for i in range(phi.shape[0]):
            for j in range(phi.shape[1]):
                if (alpha[i,j] > phi[i,j]):
                    theta[i,j] =  (alpha[i,j] - phi[i,j])
                else:
                    theta[i,j] =  2.0*np.pi + (alpha[i,j] - phi[i,j])


        magFlpara =  magpol_p*np.cos(theta)
        Flpara_r_p = magFlpara*np.cos(phi)
        Flpara_z_p = magFlpara*np.sin(phi)
        magFlpara_p = np.sqrt(Flpara_r_p*Flpara_r_p + Flpara_z_p*Flpara_z_p)
        
        f1 = plt.figure(num=1)
        ax1 = f1.add_subplot(111)
        ax1.set_aspect(1.0)
        plt.quiver(x1new,x2new,(Flpara_r_p/magFlpara_p).T,(Flpara_z_p/magFlpara_p).T,color=kwargs.get('color','k'),headwidth=kwargs.get('hw',5),width=kwargs.get('width',0.001),scale_units='xy')
        I.myfieldlines(D,[3.0,6.0,10.0,15.0],[0.01,0.01,0.01,0.01],stream=True,colors='k')

    #FOR THE PERPENDICULAR FORCES WE HAVE TO PROJECT THE Fl TO THE Bp
    if kwargs.get('Perp',False)==True:
        phi = np.arctan2(b2new,b1new)
        theta = np.zeros(phi.shape)
        Flper_r_p=np.zeros(phi.shape)
        Flper_z_p=np.zeros(phi.shape)
        magFlper=np.zeros(phi.shape)
        alpha = np.arctan2(Flz_p,Flr_p)


        for i in range(phi.shape[0]): 
            for j  in range(phi.shape[1]):
                if (alpha[i,j] > phi[i,j]):
                    theta[i,j] =  (alpha[i,j] - phi[i,j])
                    magFlper[i,j] = magpol_p[i,j]*np.sin(theta[i,j])
                    Flper_r_p[i,j] = -1.0*magFlper[i,j]*np.sin(phi[i,j])
                    Flper_z_p[i,j] = magFlper[i,j]*np.cos(phi[i,j])
                else: 
                    theta[i,j] =  2.0*np.pi + (alpha[i,j] - phi[i,j])
                    magFlper[i,j] = magpol_p[i,j]*np.sin(theta[i,j])
                    Flper_r_p[i,j] = -1.0*magFlper[i,j]*np.sin(phi[i,j])
                    Flper_z_p[i,j] = magFlper[i,j]*np.cos(phi[i,j])


        magFlper_p = np.sqrt(Flper_r_p*Flper_r_p + Flper_z_p*Flper_z_p)
        f1 = plt.figure(num=1)
        ax2 = f1.add_subplot(111)
        ax2.set_aspect(1.0)
                             
        plt.quiver(x1new,x2new,(Flper_r_p/magFlper_p).T,(Flper_z_p/magFlper_p).T,color=kwargs.get('color','k'),headwidth=kwargs.get('hw',5),width=kwargs.get('width',0.001),scale_units='xy')
        I.myfieldlines(D,[3.0,6.0,10.0,15.0],[0.01,0.01,0.01,0.01],stream=True,colors='k')


    

    
                
            
def polar2cartesian(r, t, grid, x, y, order=3):

    X, Y = np.meshgrid(x, y)

    new_r = np.sqrt(X*X+Y*Y)
    new_t = np.arctan2(X, Y)

    ir =interp1d(r, np.arange(len(r)), bounds_error=False)
    it =interp1d(t, np.arange(len(t)))

    new_ir = ir(new_r.ravel())
    new_it = it(new_t.ravel())

    new_ir[new_r.ravel() > r.max()] = len(r)-1
    new_ir[new_r.ravel() < r.min()] = 0

    return map_coordinates(grid, np.array([new_ir, new_it]),
                            order=order).reshape(new_r.shape)

##############################DISK PLOTS#############################################


def animate_image(w_dir=None,**kwargs):
    if w_dir == None: w_dir=os.getcwd()
    Vkep = 1.0e-5*np.sqrt((phc.G*phc.Msun*30.0)/(phc.au))
    frnum = kwargs.get('frames',plp.time_info(w_dir=w_dir)['nlast'])

    D0 = plp.pload(0,w_dir=w_dir)
    
    
    cdir = os.getcwd()
    files = []
    os.chdir(w_dir)
    os.system('mkdir ysymv2_movie')
   
    fig = plt.figure(num=1,figsize=[10,10])
    ax = fig.add_subplot(111)
   
    
    for i in range(frnum):
        D = plp.pload(i,w_dir=w_dir)
        ax.cla()
        
        
        if len(fig.axes) > 1:
            fig.delaxes(fig.axes[1])
            fig.subplots_adjust(right=0.90)
        
        [plc1,plc2,plcol,plq1,plq2,plfl1,plfl2,plfl3,plfl4] =symimage(D,'v2',ysym=True,flines=(True,[3.0,5.0,8.0,12.0,15.0],[0.01,0.01,0.01,0.01,0.01]),scalefact=Vkep,cbar=(True,'vertical'),arrows=(True,10),label1 = r'Radius [AU]', label2=r'Height [AU]',title=r'$\rm V_{\rm jet}$ [km s$^{-1}$]')
        ax.set_aspect('equal')
        
        axcolor = 'white'
        axtime = plt.axes([0.18, 0.01, 0.55, 0.02], axisbg=axcolor)
        axtime.cla()
        stime = Slider(axtime, 'Nsteps', 0, frnum, valinit=i)
        
        
        fname = w_dir+'ysymv2_movie/_tmp%03d.png'%i
        print 'Saving frame', fname
        fig.savefig(fname)
        files.append(fname)
        
        fig.delaxes(axtime)
        
    
    

    print 'Making movie animation.avi - this make take a while'
    os.system("mencoder 'mf://ysymv2_movie/_tmp*.png' -mf type=png:fps=15 -ovc lavc -lavcopts vcodec=msmpeg4:vbitrate=16000  -o ysymv2_movie/animation.avi")
    os.chdir(cdir)
    
        
def animate_plot(w_dir=None, **kwargs):
    if w_dir == None: w_dir=os.getcwd()

    cdir = os.getcwd()
    files = []
    os.chdir(w_dir)
    os.system('mkdir j_movie')
    D0 = plp.pload(0,w_dir=w_dir)
    fig = plt.figure(num=1,figsize=[6,12])
    ax = fig.add_subplot(111)
    frnum = kwargs.get('frames',D0.time_info()['nlast'])
    
    for i in range(frnum):  # 50 frames
        D = plp.pload(i,w_dir=w_dir)
        Ra = da.Rad_Average()
        xitem = D.x1
	yitem = D.x1*D.v3[:,32,10]*D.rho[:,32,10]
    #    yitem = Ra.Sigma(D,ul=1.0,urho=1.0e-9,Mstar=10.0,Gammae=5.0/3.0)
        ax.cla()
        
        
        

        if kwargs.get('pltype','normal') == 'normal':
            ax.plot(xitem,yitem,'k-')
        if kwargs.get('pltype','normal') == 'logx':
            ax.semilogx(xitem,yitem,'k-')
        if kwargs.get('pltype','normal') == 'logy':
            ax.semilogy(xitem,yitem,'k-')
        if kwargs.get('pltype','normal') == 'loglog':
            ax.loglog(xitem,yitem,'k-')

        ax.minorticks_on()
        ax.set_xlabel(kwargs.get('xlabel',r'XLabel'))
        ax.set_ylabel(kwargs.get('ylabel',r'YLabel'))
        ax.set_title(kwargs.get('title',r'Title'))
        
        ax.text(90.0,1.3,r'$N_{\rm rot} = %04d$'%i)
        fname = w_dir+'j_movie/_tmp%03d.png'%i
        
        print 'Saving frame', fname
        fig.savefig(fname)
        files.append(fname)
        
    

    print 'Making movie animation.mpg - this make take a while'
    os.system("mencoder 'mf://j_movie/_tmp*.png' -mf type=png:fps=10 -ovc lavc -lavcopts vcodec=wmv2 -oac copy -o movie/animation.mpg")
    os.chdir(cdir)






    





