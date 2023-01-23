#import project1 as p1
#import utils
import math
import csv
from mpl_toolkits.mplot3d import Axes3D
import powerlaw

## Other main packages
from astropy.time import Time
import astropy.units as u
## Usual packages
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
import os
from os.path import exists

from astroquery.jplhorizons import Horizons
import sys
global PYTHON3

from astropy.time import Time
from astroquery.jplhorizons import Horizons
import os
from io import StringIO

sim_start_date = "2000-01-01.5"     # simulating a solar system starting from this date
sim_duration = 2 * 365                # (int) simulation duration in days

time = Time(sim_start_date).jd
SB_Constant = 5.67*10**(-8)
L_Sol = 3.828*10**26
G = 6.672*10**(-11)
L_SOL_Earth=1.75*10**17
K_2 = 0.25
QualityF = 10

B_Constant = 1.38*10**(-23)

def rndmPLaw(a, b, g, size=1):
    """Power-law gen for pdf(x)\propto x^{g-1} for a<=x<=b"""
    r = np.random.random(size=size)
    ag, bg = a**g, b**g
    return (ag + (bg - ag)*r)**(1./g)

def plot_tempVsSemiMajorFlux(semiMajor, temperatureMoons):

    plt.subplots()

    npSemiMajor = np.array(semiMajor)
    npTemperature = np.array(temperatureMoons)
    plt.scatter(npSemiMajor, npTemperature, s=5, c='b')

    xs = npSemiMajor #semiMajor
    ys = npTemperature #temperatureMoons


    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.xlabel("Semi Major Axis [AU]")
    plt.ylabel("Equilibrium Temperature of Solar System Moons [K]")
    plt.show()

def plot_tempVsSemiMajorFluxPLaw(semiMajor, temperatureMoons, tempSolMoons, lsSolarSMMoon):

    plt.subplots()

    npSemiMajor = np.array(semiMajor)
    npTemperature = np.array(temperatureMoons)
    plt.scatter(npSemiMajor, npTemperature, s=2, c='black')
    plt.scatter(lsSolarSMMoon, tempSolMoons, s=5, c='blue')

    xs = npSemiMajor #semiMajor
    ys = npTemperature #temperatureMoons


    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    plt.xlabel("Semi Major Axis [AU]")
    plt.ylabel("Equilibrium Temperature of Solar System Moons [K]")
    plt.show()

def plot_tempVsSemiMajorTidal(semiMajor, tempertaurTidalHeating, K_2Qratio):

    plt.subplots()

    npSemiMajor = np.array(semiMajor)
    npTemperature = np.array(tempertaurTidalHeating)

    plt.scatter(npSemiMajor, npTemperature, s=5, c='b')
    xs = npSemiMajor #semiMajor
    ys = npTemperature #temperatureMoons

    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    #plt.text(0.1, 230, 'K2/Q Ratio =' + "{:.2f}".format(K_2Qratio), color='red', size='medium')
    plt.xlabel("Semi Major Axis [AU]")
    plt.ylabel("Equilibrium Temperature of Solar System Moons [K]")
    plt.show()

def plot_tempVsSemiMajor(semiMajor, temperatureMoons, tempertaurTidalHeating, K_2Qratio):

    fig, ax = plt.subplots(2)
    index = 0
    colors = []
    for var in temperatureMoons:
        index += 1
        if(index <= 3):
            colors.append('b')
        else:
            if (index > 3 and index <= 13):
                colors.append('r')
            else:
                if (index > 13 and index <= 30):
                    colors.append('g')
                else:
                    if (index > 30 and index <= 100):
                        colors.append('y')
                    else:
                        colors.append('orange')

    npSemiMajor = np.array(semiMajor)
    npTemperature = np.array(temperatureMoons)
    npTidalTemperature = np.array(tempertaurTidalHeating)
    ax[0].scatter(npSemiMajor, npTemperature, s=5, c='b')
    ax[1].scatter(npSemiMajor, npTemperature, s=5, c='b')

    xs = npSemiMajor #semiMajor
    ys = npTemperature #temperatureMoons

    scale_factor = 1.5

    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    if (len(npTidalTemperature) >0):
        ax[1].scatter(npSemiMajor, npTidalTemperature, s=5, c='b')
        xmin, xmax = ax[1].axis()[:2]
        xs = npSemiMajor  # semiMajor
        ys = npTidalTemperature  # temperatureMoons

    #ax[1].text(0.1, 230, 'K2/Q Ratio ='+"{:.2f}".format(K_2Qratio), color='red',size='medium')

    ax[0].set_xlabel("Semi Major Axis [AU]",size='medium')
    ax[0].set_ylabel("Equilibrium Temperature of Solar System Moons [K]",size='medium')
    ax[1].set_xlabel("Semi Major Axis [AU]",size='medium')
    ax[1].set_ylabel("Equilibrium Temperature of Solar System Moons [K]",size='medium')
    plt.show()

def plot_tempStarVsSemiMajor(semiMajor, temperature):

    plt.subplots()
    index = 0
    colors = []
    npSemiMajor = np.array(semiMajor)
    npTemperature = np.array(temperature)

    colors = []
    index = 0
    for var in npTemperature:
        if (npTemperature[index] == 5778):  # Sun
            colors.append("red")
        else:
            colors.append("black")
        index = index + 1
    plt.scatter(npSemiMajor, npTemperature, s=5, c = colors)
    xmin, xmax = plt.axis()[:2]
    xs = npSemiMajor #semiMajor
    ys = npTemperature #temperatureMoons
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.xlabel("Semi Major Axis of Moons [AU]")
    plt.ylabel("Temperature of Host Star [K]")
    plt.show()

def plot_tempStarVsSemiMajorPlanet(semiMajor, temperatureBig, temperatureSmall, lbHabitable, ubHabitable,lstSolarPlanet):

    plt.subplots()
    index = 0
    colors = []
    npSemiMajor = np.array(semiMajor)
    npTemperature = np.array(temperatureBig)
    nptemperatureSmall = np.array(temperatureSmall)

    colors = []
    index = 0
    for var in npTemperature:
        if (npTemperature[index] == 5778):  # Sun
            colors.append("red")
        else:
            colors.append("black")
        index = index + 1

    habitableLB = lbHabitable
    habitableUB = ubHabitable
    starTemperature = npTemperature
    lstSun = []
    for lstSM in lstSolarPlanet:
        lstSun.append(5778)
    plt.scatter(npSemiMajor, npTemperature, s=2, c=colors)
    plt.scatter(lstSolarPlanet, lstSun, s=10, c='red')
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()
    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    nptemperatureSmall.sort()
    habitableLB.sort()
    habitableUB.sort()

    index = 0
    for temp in nptemperatureSmall:
        plt.axhline(temp, xmin = 0.035 + habitableLB[index]/((xmax-xmin)), xmax =  0.035 + habitableUB[index]/((xmax-xmin)), color='lightgrey', markersize=1)
        index+=1
    plt.scatter(lstSolarPlanet, lstSun, s=10, c='red')
    plt.xlabel("Semi Major Axis of Planets [AU]")
    plt.ylabel("Temperature of Host Star [K]")
    plt.show()

def plot_tempStarVsSemiMajorPlanetPLaw(semiMajor, temperatureBig, temperatureSmall, lbHabitable, ubHabitable,lstSolarPlanet):

    plt.subplots()
    index = 0
    colors = []
    npSemiMajor = np.array(semiMajor)
    npTemperature = np.array(temperatureBig)
    nptemperatureSmall = np.array(temperatureSmall)

    semiMajorAct=[]
    temperatureBigAct=[]
    temperatureSmallAct=[]
    index = 0
    for i in range(0,50):
        for j in range (0,40):
            semiMajorAct.append(semiMajor[index])
            temperatureBigAct.append(temperatureBig[index])
            temperatureSmallAct.append(temperatureSmall[index])
            index+=1
    npSemiMajorAct = np.array(semiMajorAct)
    npTemperatureAct = np.array(temperatureBigAct)
    nptemperatureSmallAct = np.array(temperatureSmallAct)

    colors = []
    index = 0
    for var in npTemperatureAct:
        if (npTemperatureAct[index] == 5778):  # Sun
            colors.append("red")
        else:
            colors.append("black")
        index = index + 1

    habitableLB = lbHabitable
    habitableUB = ubHabitable
    starTemperature = npTemperatureAct
    lstSun = []
    for lstSM in lstSolarPlanet:
        lstSun.append(5778)
    plt.scatter(npSemiMajorAct, npTemperatureAct, s=2, c=colors)
    plt.scatter(lstSolarPlanet, lstSun, s=10, c='red')
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)

    index = 0
    for temp in nptemperatureSmall:
        plt.axhline(temp, xmin = 0.035 + habitableLB[index]/((xmax-xmin)), xmax =  0.035 + habitableUB[index]/((xmax-xmin)), color='lightgrey', markersize = 1)
        index+=1
    plt.scatter(lstSolarPlanet, lstSun, s=10, c='red')
    plt.xlabel("Semi Major Axis of Planets [AU]")
    plt.ylabel("Temperature of Host Star [K]")
    plt.show()

def plot_tempStarVsSemiMajorMoonPLaw(semiMajor, temperatureBig, temperatureSmall,lstSolarMoon):

    plt.subplots()
    index = 0
    colors = []
    npSemiMajor = np.array(semiMajor)
    npTemperature = np.array(temperatureBig)
    nptemperatureSmall = np.array(temperatureSmall)

    colors = []
    index = 0
    for var in npTemperature:
        if (npTemperature[index] == 5778):  # Sun
            colors.append("red")
        else:
            colors.append("black")
        index = index + 1

    starTemperature = npTemperature
    lstSun = []
    for lstSM in lstSolarMoon:
        lstSun.append(5778)
    plt.scatter(npSemiMajor, npTemperature, s=2, c=colors)
    plt.scatter(lstSolarMoon, lstSun, s=5, c='blue')
    xmin, xmax = plt.xlim()
    ymin, ymax = plt.ylim()

    plt.xlim(xmin, xmax)
    plt.ylim(ymin, ymax)
    plt.scatter(lstSolarMoon, lstSun, s=10, c='blue')
    plt.xlabel("Semi Major Axis of Moons [AU]")
    plt.ylabel("Temperature of Host Star [K]")
    plt.show()

def powerFromAlbedoL(distance, radius, albedo):
    solar_constant = TSI(distance)
    albedoconv = float(albedo)
    rConv = float(radius)
    power = solar_constant*(1-albedoconv)*math.pi*(rConv**2)
    return power

def powerFromAlbedo(solar_constant, albedo):
    albedoconv = float(albedo)
    power = ((solar_constant*(1-albedoconv)/(4*SB_Constant)))
    return power

def tempFromAlbedo(solar_constant, albedo):
    albedoconv = float(albedo)
    temperature = ((solar_constant*(1-albedoconv)/(4*SB_Constant)))**0.25
    return temperature

def TSI(R):
    rConv = float(R)*1.496*(10)**8*1000
    TSI = (L_Sol)/(4*(math.pi)*rConv**(2))
    return TSI

def TSI1(R):
    rConv = R*1.496*(10)**8*1000
    TSI = (L_Sol)/(4*(math.pi)*rConv**(2))
    return TSI

def tempFromAlbedo1(distance, albedo):
    rConv = float(distance)
    albedoconv = float(albedo)
    return (280*((1-albedoconv)/rConv)**0.25)

def tempOfStar(radius):
    solar_constant = TSI(radius)
    print("TSI",solar_constant)
    temperature = ((solar_constant/(SB_Constant)))**0.25
    return temperature

def isfloat(num):
    try:
        float(num)
        return True
    except ValueError:
        return False

def findNonSpace(strText,iIndexEq):
    for i in range(iIndexEq, len(strText)-1):
        if strText[i] != " " :
            iIndexStart = i
            iIndexEnd = strText.find(" ", iIndexStart+1,len(strText)-1)
            valRadius = strText[iIndexStart:iIndexEnd-1]
            return valRadius
    valRadius = "50"
    return valRadius

def find_MeanRadius(strText):
    iIndex = strText.find("radius")
    if (iIndex == -1):
        iIndex = strText.find("Radius")
    if (iIndex ==-1):
        return "50",False
    else:
        iStart = iIndex
        iIndexEq = strText.find("=",iStart,len(strText)-1)
        if(iIndexEq == -1):
            valRadius = "50"
            return "50", False
        else:
            iIndex = strText.find("+", iIndexEq, len(strText) - 1)
            if iIndex > -1:
                if ( iIndex - iIndexEq >= 12):
                    iIndex = -1
            if (iIndex == -1):
                iIndex = strText.find(".", iIndexEq, len(strText) - 1)
                if ( iIndex - iIndexEq >= 12):
                    iIndex = -1
                if (iIndex == -1):
                    iIndex = strText.find("x", iIndexEq, len(strText) - 1)
                    if (iIndex == -1):
                        iIndex = strText.find("Density", iIndexEq, len(strText) - 1)
                        if (iIndex == -1):
                            valRadius = "50"
                            return "50", False
                        else:
                            valRadius = strText[iIndexEq + 1:iIndex]
                            valRadius = valRadius.strip()
                            return valRadius, True
                    else:
                        valRadius = strText[iIndexEq + 1:iIndex]
                        valRadius = valRadius.strip()
                        return valRadius, True
                else:
                    valRadius = strText[iIndexEq + 1:iIndex+2]
                    valRadius=valRadius.strip()
                    return valRadius, True
            else:
                valRadius = strText[iIndexEq+1:iIndex-1]
                valRadius = valRadius.strip()
                return valRadius, True
    return valRadius,True

def find_albedo(strText):
    """
    Loads the 2D toy dataset as numpy arrays.
    Returns the tuple (features, labels) in which features is an Nx2 numpy matrix and
    labels is a length-N vector of +1/-1 labels.
    """
    iIndex = strText.find("Albedo")
    if (iIndex ==-1):
        return "0.04",False

    else:
        iStart = iIndex
        iIndex = strText.find("=",iStart,len(strText)-1)
        if(iIndex == -1):
            valAlbedo = "0.04"
            return "0.04", False
        else:
            iIndex = strText.find(".", iStart, len(strText) - 1)
            if (iIndex == -1):
                valAlbedo = "0.04"
                return "0.04", False
            else:
                valAlbedo = strText[iIndex-1:iIndex+3]
                return valAlbedo, True
    return valAlbedo,True

def find_SatelliteFileName(strText):
    """
    Loads the 2D toy dataset as numpy arrays.
    Returns the tuple (features, labels) in which features is an Nx2 numpy matrix and
    labels is a length-N vector of +1/-1 labels.
    """
    iIndex = strText.find("(")
    if (iIndex ==-1):
        return "", False
    else:
        strSatellite = strText[0:iIndex-1]
        return strSatellite,True

def load_dataReturnMeanRad(strFileName):
    path_data = "..\Habitability\\"
    path_data1 = ".txt"
    path_data = path_data + strFileName + path_data1

    file_exists = exists(path_data)

    if sys.version_info[0] < 3:
        PYTHON3 = False
    else:
        PYTHON3 = True

    if PYTHON3:
        if(file_exists):
            f_data = open(path_data, encoding="latin1")
        else:
            valRadius = "50"
            return valRadius
    else:
        if(file_exists):
            f_data = open(path_data)
        else:
            valRadius = "valRadius"
            return valRadius

    strListContents = f_data.readlines()
    bRadiusFound = False
    for content in strListContents:
        valRadius, bFound = find_MeanRadius(content)
        if (bFound == True):
            bRadiusFound = True
            valRadius = valRadius.strip()
            valRadius.replace(" ","")
            print("Radius before =", valRadius)
            if(isfloat(valRadius)):
                return valRadius
            else:
                valRadius = "50"
                return valRadius
    if(bRadiusFound== False):
       valRadius = "50"
    f_data.close()
    return valRadius

def load_dataReturnAlbedo(strFileName):
    path_data = "..\Habitability\\"
    path_data1 = ".txt"
    path_data = path_data + strFileName + path_data1

    file_exists = exists(path_data)

    if sys.version_info[0] < 3:
        PYTHON3 = False
    else:
        PYTHON3 = True

    if PYTHON3:
        if(file_exists):
            f_data = open(path_data, encoding="latin1")
        else:
            valAlbedo = "0.04"
            return valAlbedo
    else:
        if(file_exists):
            f_data = open(path_data)
        else:
            valAlbedo = "0.04"
            return valAlbedo

    strListContents = f_data.readlines()
    bAlbedoFound = False
    for content in strListContents:
        valAlbedo, bFound = find_albedo(content)
        if (bFound == True):
            bAlbedoFound = True
            valAlbedo = valAlbedo.strip()
            valAlbedo.replace(" ","")
            if(isfloat(valAlbedo)):
                return valAlbedo
            else:
                valAlbedo = "0.04"
                return valAlbedo
    if(bAlbedoFound== False):
       valAlbedo = "0.04"
    f_data.close()
    return valAlbedo

def BuildSatelliteDB():
    lstMoon = []
    lstEphemeris = []
    lstAlbedo = []
    lstMeanRad = []
    lstSemiMajor = []
    lstDistHostStar = []
    lstTemp = []
    lstPlanet = []
    lstMPlanetT = []
    lstEccentricty = []
    lstMPlanet = [5.9736*(10)**(24),6.39*(10)**(23),1.898*(10)**(27),5.683*(10)**(26),8.681*(10)**(25),1.024*(10)**(26)]
    for i in range(3,9,1):
        lstMoonData = []
        lstEph = []
        lstVector = []
        lstAlbedoSat = []
        lstMeanRadSat = []
        lstSemiMajorObj = []
        lstTempObj = []
        lstDistStarObj = []
        lstMPlanetObj = []
        lstEccentricityObj=[]
        if (i==3): #Earth
          iIndex = i*100+1
          objSemiMajor = Horizons(id=399, location='@Sun', epochs=time, id_type='id')
          lstPlanet.append(objSemiMajor.elements()[0]['a'])
          for iMoon in range(iIndex, iIndex+1, 1):
              obj = Horizons(id=iMoon, location=399, epochs=time, id_type='id')
              el = obj.elements()
              distSemiMajor = el[0]['a']
              eccentricity = el[0]['e']
              distance= objSemiMajor.elements()[0]['a']+ el[0]['a']
              lstDistStarObj.append(distance)
              lstSemiMajorObj.append(distSemiMajor)
              strStatelliteName,bFound = find_SatelliteFileName(el[0][0])
              print("Corrected Satellite = ",strStatelliteName)
              valAlbedo = load_dataReturnAlbedo(strStatelliteName)
              print("Albedo Main= ",valAlbedo)
              valRadius = load_dataReturnMeanRad(strStatelliteName)
              tempCalculated = tempFromAlbedo(TSI(distance),valAlbedo)
              #tempCalculated = tempFromAlbedo1(distance, valAlbedo)
              #tempCalculated = tempFromAlbedo1(distSemiMajor, valAlbedo)
              eph = obj.ephemerides()
              lstMoonData.append(el)
              lstEph.append(eph)
              lstAlbedoSat.append(valAlbedo)
              lstMeanRadSat.append(valRadius)
              lstTempObj.append(tempCalculated)
              lstMPlanetObj.append(lstMPlanet[0])
              lstEccentricityObj.append(eccentricity)
        if (i==4): #Mars
          iIndex = i*100+1
          objSemiMajor = Horizons(id=499, location='@Sun', epochs=time, id_type='id')
          lstPlanet.append(objSemiMajor.elements()[0]['a'])
          for iMoon in range(iIndex, iIndex+2, 1):
              obj = Horizons(id=iMoon, location=499, epochs=time , id_type='id')
              epochs={'start':'2017-01-01',
                         'stop':'2017-02-01',
                         'step':'1d'}
              el = obj.elements()
              distance= objSemiMajor.elements()[0]['a']+ el[0]['a']
              distSemiMajor = el[0]['a']
              lstDistStarObj.append(distance)
              lstSemiMajorObj.append(distSemiMajor)
              strStatelliteName, bFound = find_SatelliteFileName(el[0][0])
              eccentricity = el[0]['e']
              print("Corrected Satellite = ", strStatelliteName)
              valAlbedo = load_dataReturnAlbedo(strStatelliteName)
              print("Albedo Main= ",valAlbedo)
              valRadius = load_dataReturnMeanRad(strStatelliteName)
              tempCalculated = tempFromAlbedo(TSI(distance),valAlbedo)
              #tempCalculated = tempFromAlbedo1(distance, valAlbedo)
              #tempCalculated = tempFromAlbedo1(distSemiMajor, valAlbedo)
              eph = obj.ephemerides()
              lstMoonData.append(el)
              lstEph.append(eph)
              lstAlbedoSat.append(valAlbedo)
              lstMeanRadSat.append(valRadius)
              lstTempObj.append(tempCalculated)
              lstMPlanetObj.append(lstMPlanet[1])
              lstEccentricityObj.append(eccentricity)
        if (i== 5): #'Jupiter':
          iIndex = i*100+1
          objSemiMajor = Horizons(id=599, location='@Sun', epochs=time, id_type='id')
          lstPlanet.append(objSemiMajor.elements()[0]['a'])
          for iMoon in range(iIndex, iIndex+71, 1):
              obj = Horizons(id=iMoon, location=599, epochs=time , id_type='id')
              el = obj.elements()
              distance= objSemiMajor.elements()[0]['a']+ el[0]['a']
              distSemiMajor = el[0]['a']
              lstDistStarObj.append(distance)
              lstSemiMajorObj.append(distSemiMajor)
              strStatelliteName, bFound = find_SatelliteFileName(el[0][0])
              eccentricity = el[0]['e']
              print("Corrected Satellite = ", strStatelliteName)
              valAlbedo = load_dataReturnAlbedo(strStatelliteName)
              print("Albedo Main= ",valAlbedo)
              valRadius = load_dataReturnMeanRad(strStatelliteName)
              tempCalculated = tempFromAlbedo(TSI(distance),valAlbedo)
              #tempCalculated = tempFromAlbedo1(distance, valAlbedo)
              #tempCalculated = tempFromAlbedo1(distSemiMajor, valAlbedo)
              eph = obj.ephemerides()
              lstMoonData.append(el)
              lstEph.append(eph)
              lstAlbedoSat.append(valAlbedo)
              lstMeanRadSat.append(valRadius)
              lstTempObj.append(tempCalculated)
              lstMPlanetObj.append(lstMPlanet[2])
              lstEccentricityObj.append(eccentricity)
        if (i== 6): #'Saturn'
          iIndex = i*100+1
          objSemiMajor = Horizons(id=699, location='@Sun', epochs=time, id_type='id')
          lstPlanet.append(objSemiMajor.elements()[0]['a'])
          for iMoon in range(iIndex, iIndex+66, 1):
              obj = Horizons(id=iMoon, location=699, epochs=time, id_type='id')
              el = obj.elements()
              distance= objSemiMajor.elements()[0]['a']+ el[0]['a']
              distSemiMajor = el[0]['a']
              lstDistStarObj.append(distance)
              lstSemiMajorObj.append(distSemiMajor)
              strStatelliteName, bFound = find_SatelliteFileName(el[0][0])
              eccentricity = el[0]['e']
              print("Corrected Satellite = ", strStatelliteName)
              valAlbedo = load_dataReturnAlbedo(strStatelliteName)
              print("Albedo Main= ",valAlbedo)
              valRadius = load_dataReturnMeanRad(strStatelliteName)
              tempCalculated = tempFromAlbedo(TSI(distance),valAlbedo)
              #tempCalculated = tempFromAlbedo1(distance, valAlbedo)
              #tempCalculated = tempFromAlbedo(distance, valAlbedo)
              eph = obj.ephemerides()
              lstMoonData.append(el)
              lstEph.append(eph)
              lstAlbedoSat.append(valAlbedo)
              lstMeanRadSat.append(valRadius)
              lstTempObj.append(tempCalculated)
              lstMPlanetObj.append(lstMPlanet[3])
              lstEccentricityObj.append(eccentricity)
        if (i== 7): #'Uranus'
          iIndex = i*100+1
          objSemiMajor = Horizons(id=799, location='@Sun', epochs=time, id_type='id')
          lstPlanet.append(objSemiMajor.elements()[0]['a'])
          for iMoon in range(iIndex, iIndex+24, 1):
              obj = Horizons(id=iMoon, location=799, epochs=time, id_type='id')
              el = obj.elements()
              distance= objSemiMajor.elements()[0]['a']+ el[0]['a']
              distSemiMajor = el[0]['a']
              lstDistStarObj.append(distance)
              lstSemiMajorObj.append(distSemiMajor)
              strStatelliteName, bFound = find_SatelliteFileName(el[0][0])
              eccentricity = el[0]['e']
              print("Corrected Satellite = ", strStatelliteName)
              valAlbedo = load_dataReturnAlbedo(strStatelliteName)
              print("Albedo Main= ",valAlbedo)
              valRadius = load_dataReturnMeanRad(strStatelliteName)
              tempCalculated = tempFromAlbedo(TSI(distance),valAlbedo)
              #tempCalculated = tempFromAlbedo1(distance, valAlbedo)
              #tempCalculated = tempFromAlbedo1(distSemiMajor, valAlbedo)
              eph = obj.ephemerides()
              lstMoonData.append(el)
              lstEph.append(eph)
              lstAlbedoSat.append(valAlbedo)
              lstMeanRadSat.append(valRadius)
              lstTempObj.append(tempCalculated)
              lstMPlanetObj.append(lstMPlanet[4])
              lstEccentricityObj.append(eccentricity)
        if (i== 8): #'Neptune'
          iIndex = i*100+1
          objSemiMajor = Horizons(id=899, location='@Sun', epochs=time, id_type='id')
          lstPlanet.append(objSemiMajor.elements()[0]['a'])
          for iMoon in range(iIndex, iIndex+14, 1):
              obj = Horizons(id=iMoon, location=899, epochs=time, id_type='id')
              el = obj.elements()
              distance= objSemiMajor.elements()[0]['a']+ el[0]['a']
              distSemiMajor = el[0]['a']
              lstDistStarObj.append(distance)
              lstSemiMajorObj.append(distSemiMajor)
              strStatelliteName, bFound = find_SatelliteFileName(el[0][0])
              eccentricity = el[0]['e']
              print("Corrected Satellite = ", strStatelliteName)
              valAlbedo = load_dataReturnAlbedo(strStatelliteName)
              print("Albedo Main= ",valAlbedo)
              valRadius = load_dataReturnMeanRad(strStatelliteName)
              tempCalculated = tempFromAlbedo(TSI(distance),valAlbedo)
              #tempCalculated = tempFromAlbedo1(distance, valAlbedo)
              #tempCalculated = tempFromAlbedo1(distSemiMajor, valAlbedo)
              eph = obj.ephemerides()
              lstMoonData.append(el)
              lstEph.append(eph)
              lstAlbedoSat.append(valAlbedo)
              lstMeanRadSat.append(valRadius)
              lstTempObj.append(tempCalculated)
              lstMPlanetObj.append(lstMPlanet[5])
              lstEccentricityObj.append(eccentricity)
        if (i== 9): #'Pluto'
          iIndex = i*100+1
          objSemiMajor = Horizons(id=999, location='@Sun', epochs=time, id_type='id')
          lstPlanet.append(objSemiMajor.elements()[0]['a'])
          for iMoon in range(iIndex, iIndex+6, 1):
              obj = Horizons(id=iMoon, location=999, epochs=time, id_type='id')
              el = obj.elements()
              distance= float(objSemiMajor.elements()[0]['a']+ el[0]['a'])
              distSemiMajor = el[0]['a']
              lstDistStarObj.append(distance)
              lstSemiMajorObj.append(distSemiMajor)
              strStatelliteName, bFound = find_SatelliteFileName(el[0][0])
              eccentricity = el[0]['e']
              print("Corrected Satellite = ", strStatelliteName)
              valAlbedo = float(load_dataReturnAlbedo(strStatelliteName))
              print("Albedo Main= ",valAlbedo)
              valRadius = float(load_dataReturnMeanRad(strStatelliteName))
              tempCalculated = tempFromAlbedo(TSI(distance),valAlbedo)
              #tempCalculated = tempFromAlbedo1(distance, valAlbedo)
              #tempCalculated = tempFromAlbedo1(distSemiMajor, valAlbedo)
              eph = obj.ephemerides()
              lstMoonData.append(el)
              lstEph.append(eph)
              lstAlbedoSat.append(valAlbedo)
              lstMeanRadSat.append(valRadius)
              lstTempObj.append(tempCalculated)
              lstMPlanetObj.append(lstMPlanet[6])
              lstEccentricityObj.append(eccentricity)
        lstMoon.append(lstMoonData)
        lstEphemeris.append(lstEph)
        lstAlbedo.append(lstAlbedoSat)
        lstMeanRad.append(lstMeanRadSat)
        lstSemiMajor.append(lstSemiMajorObj)
        lstDistHostStar.append(lstDistStarObj)
        lstTemp.append(lstTempObj)
        lstMPlanetT.append(lstMPlanetObj)
        lstEccentricty.append(lstEccentricityObj)
    return lstMoon, lstEphemeris, lstAlbedo,lstMeanRad,lstSemiMajor,lstTemp,lstDistHostStar,lstPlanet,lstMPlanetT,lstEccentricty

lstMoon, lstEphemeris, lstAlbedo,lstMeanRad,lstSemiMajor,lstTemp,lstDistHostStar,lstPlanet,lstMPlanetT,lstEccentricty = BuildSatelliteDB()

# Examples How to retrieve data from Ephemeris V Brightness illumination
print (' Stored Radius',lstMeanRad)
tempSt = tempOfStar(0.00465047) #(" Temperature of sun",tempSt )


import random

# Random float number
lstSM1D = []

lstTemp1D = []
lstHostStarTemp = []
lstDistStar1D = []
lstPlanet1D = []
lstHostStarPlanet = []
lstMeanRad1D = []
lstTempStar = []
lstEccentricty1D = []
lstMPlanetT1D = []
lstAlb1d = []

# Randomly generated Host Star Temp (We need to apply Power Law)

for i in range(50):
    temp = random.uniform(2000, 20000)
    lstTempStar.append(temp)

lstTempStar.append(5778)
lstTempStar.sort()

for lstSM in lstSemiMajor:
    for value in lstSM:
        lstSM1D.append(value)
        lstHostStarTemp.append(tempSt)

for lsttempAlb in lstAlbedo:
    for value in lsttempAlb:
        lstAlb1d.append(value)

for lstTempVal in lstTemp:
    for value in lstTempVal:
        lstTemp1D.append(value)

for lstTempMeanRad in lstMeanRad:
    for value in lstTempMeanRad:
        lstMeanRad1D.append(value)

for lstTempEccentricity in lstEccentricty:
    for value in lstTempEccentricity:
        lstEccentricty1D.append(value)

for lstTempMPlanetT in lstMPlanetT:
    for value in lstTempMPlanetT:
        lstMPlanetT1D.append(value)


for lstDist in lstDistHostStar:
    for value in lstDist:
        lstDistStar1D.append(value)

print("Planet T 1D= ", lstMPlanetT1D)


lstTidalTemp = []

index = 0

# Stellar Flux and Tidal Heating are added
for tempSM in lstSM1D:
    mplanet = lstMPlanetT1D[index]
    print("Mplanet ", index,"=", mplanet)
    radius = float(lstMeanRad1D[index]) * 1000
    eccentricity = float(lstEccentricty1D[index])
    semimajor = float(tempSM) * 1.496 * ((10) ** (8)) * 1000
    K_2 = random.uniform(0.4, 0.6)
    QualityF = random.uniform(0.005, 0.006)
    lTidal = ((21 / 2) * (K_2) * ((G) ** (3 / 2)) * ((mplanet) ** (5 / 2)) * ((radius) ** (5)) *
              ((eccentricity) ** (2))) / (QualityF * ((semimajor) ** (15 / 2)))
    pTidal = (lTidal / (4 * math.pi * SB_Constant * ((radius) ** (2))))
    distance = float(lstDistStar1D[index])
    albedo = lstAlb1d[index]
    power = powerFromAlbedo(TSI(distance),albedo)
    TTidalPlusAlbedo = (pTidal + power)**(0.25)
    lstTidalTemp.append(TTidalPlusAlbedo)
    index += 1

plot_tempVsSemiMajorFlux(lstSM1D, lstTemp1D)
plot_tempVsSemiMajorTidal(lstSM1D,lstTidalTemp, K_2/QualityF)
plot_tempVsSemiMajor(lstSM1D, lstTemp1D,lstTidalTemp, K_2/QualityF)

# Randomly generated Semi Major of Planets based on Uniform Distribution
# This is replaced by Power Law

lstSolarPlanet = lstPlanet

for i in range(50):
    lstPlanet.append(random.uniform(0.2, 30.5))

index = 0
for hoststartemp in lstTempStar:
    index = 0
    for distancePlanet in lstPlanet:
        if(index >= 6):
            # Randomly generated Semi Major of Planets (We need to apply Power Law)
            distancePlanet = random.uniform(0.2, 30.5)
            lstPlanet1D.append(distancePlanet)
        else:
            lstPlanet1D.append(distancePlanet)
        lstHostStarPlanet.append(hoststartemp)
        index = index + 1

iIndex = 0
for hoststartemp in lstTempStar:
    for lstSM in lstSemiMajor:
        for value in lstSM:
            iIndex+=1
            lstSM1D.append(value)
            lstHostStarTemp.append(hoststartemp)

# Plot Randomly Distributed Stars and Planets with Habitable Zones

lstHostStarPlanetnew = []
lstPlanet1dnew = []
for i in range(100):
     for j in range(70):
         temp = rndmPLaw(2500, 20000, g=-0.5, size=1)
         lstHostStarPlanetnew.append(temp[0])
         distancePlanet = rndmPLaw(0.2, 30.5, g=-0.5, size=1)
         lstPlanet1dnew.append(distancePlanet[0])

normed_hoststartemp = np.array(lstHostStarPlanetnew)
normed_semimajor = np.array(lstPlanet1dnew)
index = 0
for i in range(6):
    np.append(normed_hoststartemp,[5778])
    print ("Solar Distance = ", lstSolarPlanet[index])
    np.append(normed_semimajor,[lstSolarPlanet[index]])
    index +=1

lstHabitableUB = []
lstHabitableLB = []

for hostStarTemp in normed_hoststartemp:
    L_Star = L_Sol*(hostStarTemp/5778)**4
    habitableUBound = ((L_Star * (1 - 0.3) / (16 * math.pi * SB_Constant * (200) ** 4)) ** 0.5) / (1.496 * (10) ** 8 * 1000)
    habitableLBound = ((L_Star * (1 - 0.3) / (16 * math.pi * SB_Constant * (300) ** 4)) ** 0.5) / (1.496 * (10) ** 8 * 1000)
    lstHabitableLB.append(habitableLBound)
    lstHabitableUB.append(habitableUBound)

lsSolarPlanetNew = []
index = 0
for i in range(6):
    lsSolarPlanetNew.append(lstSolarPlanet[index])
    index +=1

plot_tempStarVsSemiMajorPlanetPLaw(normed_semimajor,normed_hoststartemp, normed_hoststartemp, lstHabitableLB, lstHabitableUB, lsSolarPlanetNew) #in AU

# Plot Randomly Distributed Stars and Moons with Habitable Zones

lstHostStarPlanetnew = []
lstMoon1dnew = []
lstTempMoon = []
for i in range(100):
     for j in range(100):
         temp = rndmPLaw(50, 20000, g=-0.5, size=1)
         lstHostStarPlanetnew.append(temp[0])
         distanceMoon = rndmPLaw(0.0005, 0.4, g=-0.5, size=1)
         lstMoon1dnew.append(distanceMoon[0])
         tempMoon = rndmPLaw(10, 500, g=-0.5, size=1) #random.uniform(0, 400)
         lstTempMoon.append(tempMoon)

normed_hoststartemp = np.array(lstHostStarPlanetnew)
normed_semimajorMoon = np.array(lstMoon1dnew)
normed_MoonTempStelFlux = np.array(lstTempMoon)

index = 0
for i in range(178):
    np.append(normed_hoststartemp,[5778])
    np.append(normed_semimajorMoon,[lstSM1D[index]])
    np.append(normed_MoonTempStelFlux,[lstTemp1D[index]])
    index +=1

lsSolarMoonNew = []
index = 0
for i in range(178):
    lsSolarMoonNew.append(lstSM1D[index])
    index +=1

# Temp Vs Semi Major of Moons
plot_tempStarVsSemiMajorMoonPLaw(normed_semimajorMoon,normed_hoststartemp, normed_hoststartemp,lsSolarMoonNew) #in AU

# Temp of moons Vs Semi Major of Moons stellar irradiance
# This one to be removed
#plot_tempVsSemiMajorFluxPLaw(normed_semimajorMoon, normed_MoonTempStelFlux, lstTemp1D, lsSolarMoonNew)

# Temp of moons Vs Semi Major of Moons stellar irradiance plus tidal heating

lstHostStarPlanetnew = []
lstMoon1dnew = []
lstTempMoon = []
for i in range(100):
     for j in range(100):
         temp = rndmPLaw(50, 20000, g=-0.5, size=1)
         lstHostStarPlanetnew.append(temp[0])
         distanceMoon = rndmPLaw(0.0005, 0.4, g=-0.5, size=1)
         lstMoon1dnew.append(distanceMoon[0])
         tempMoon = rndmPLaw(10, 700, g=-0.5, size=1) #random.uniform(0, 400)
         lstTempMoon.append(tempMoon)

normed_hoststartemp = np.array(lstHostStarPlanetnew)
normed_semimajorMoon = np.array(lstMoon1dnew)
normed_MoonTempStelFlux = np.array(lstTempMoon)

index = 0
for i in range(178):
    np.append(normed_hoststartemp,[5778])
    np.append(normed_semimajorMoon,[lstSM1D[index]])
    np.append(normed_MoonTempStelFlux,[lstTemp1D[index]])
    index +=1

lsSolarMoonNew = []
index = 0
for i in range(178):
    lsSolarMoonNew.append(lstSM1D[index])
    index +=1
#This one to be removed
#plot_tempVsSemiMajorFluxPLaw(normed_semimajorMoon, normed_MoonTempStelFlux, lstTidalTemp, lsSolarMoonNew)

lstHostStarPlanetnew = []
lstPlanet1dnew = []
normed_MoonPowStelFlux=[]
normed_MoonPowStelFluxPlusTidal = []
normed_MoonPowStelFluxMinusTidal = []
normed_MoonPowStelFluxPlusTidal2D = []
normed_MoonPowStelFluxMinusTidal2D = []
normed_MoonPowStelFlux2D = []
normed_MoonPowTidal2D = []

for i in range(100):
     for j in range(100):
         temp = rndmPLaw(2500, 20000, g=-0.5, size=1)
         lstHostStarPlanetnew.append(temp[0])
         distancePlanet = rndmPLaw(0.2, 30.5, g=-0.5, size=1)
         lstPlanet1dnew.append(distancePlanet[0])


normed_hoststartemp = np.array(lstHostStarPlanetnew)
normed_semimajor = np.array(lstPlanet1dnew)
index = 0
for i in range(6):
    np.append(normed_hoststartemp,[5778])
    print ("Solar Distance = ", lstSolarPlanet[index])
    np.append(normed_semimajor,[lstSolarPlanet[index]])
    index +=1

lstHostStarPlanetnew = []
lstMoon1dnew = []
lstTempMoon = []
index = 0

lstsemimajorplanet=[]

for i in range(100):
    distancePlanet = rndmPLaw(0.2, 30.5, g=-0.5, size=1)
    lstsemimajorplanet.append(distancePlanet[0])

lstsemimajorplanet.sort()

npSemiMajorPlanet = np.array(lstsemimajorplanet)

lstMoon2D = []
for i in range(100):
    distancePlanet = lstsemimajorplanet[i]
    lstSemiMajorMoon=[]
    for j in range(100):
         distanceMoon = rndmPLaw(0.002, 0.4, g=-0.5, size=1)
         lstSemiMajorMoon.append(distanceMoon[0])
    lstSemiMajorMoon.sort()
    lstMoon2D.append(lstSemiMajorMoon)

npSemiMajorMoon = np.array(lstMoon2D)

print("Semi Major Planet = ",npSemiMajorPlanet)
print("Semi Major Moon = ",npSemiMajorMoon)

#mplanet = random.uniform(2 * (10 ** 27), 3 * (10 ** 30))  # Close to Jupiter or higher distribution
#mplanet = random.uniform(2 * (10 ** 27), 2 * (10 ** 27))  # Close to Jupiter or higher distribution

mplanet = random.uniform(5.9 * (10 ** 24), 5.9 * (10 ** 24))  # Close to Earth's mass distribution

#radiusMoon = random.uniform(1800 * 1000, 4000 * 1000)  # Close to Jupiter or higher distribution
radiusMoon = random.uniform(1800 * 1000, 1800 * 1000)  # Close to Jupiter or higher distribution
eccentricity = random.uniform(0.003, 0.3)
K_2 = random.uniform(0.1, 0.5)b
QualityF = random.uniform(5, 10)
for i in range(100):
     normed_MoonPowStelFluxPlusTidal1D = []
     normed_MoonPowStelFluxMinusTidal1D = []
     normed_MoonPowStelFlux1D = []
     normed_MoonPowStelFluxTidal1D = []
     distancePlanet = npSemiMajorPlanet[i] #rndmPLaw(0.2, 30.5, g=-0.5, size=1)
     #mplanet = random.uniform(2 * (10 ** 27), 3 * (10 ** 30))  # Close to Jupiter or higher distribution
     #radiusMoon = random.uniform(20 * 1000, 4000 * 1000)  # Close to Jupiter or higher distribution
     #eccentricity = random.uniform(0.003, 0.3)
     #K_2 = random.uniform(0.1, 0.5)
     #QualityF = random.uniform(5, 10)
     for j in range(100):
         distanceMoon = npSemiMajorMoon[i][j]  #rndmPLaw(0.005, 0.4, g=-0.5, size=1)
         lstMoon1dnew.append(distanceMoon)
         tempMoon = rndmPLaw(10, 700, g=-0.5, size=1) #random.uniform(0, 400)
         distanceStarMoon = distancePlanet + distanceMoon
         #radiusMoon = random.uniform(20 * 1000, 4000 * 1000)  # Close to Jupiter or higher distribution
         semiMajorau = distanceMoon #random.uniform(0.1, 0.4)
         semimajor = semiMajorau * 1.496 * ((10) ** (8)) * 1000 #Jupiter Europa distribution
         #eccentricity = random.uniform(0.003, 0.3)
         power = powerFromAlbedoL(distanceStarMoon,radiusMoon, 0.3)
         normed_MoonPowStelFlux.append(power)
         #K_2 = random.uniform(0.1, 0.5)
         #QualityF = random.uniform(5, 10)
         lTidal = ((21 / 2) * (K_2) * ((G) ** (3 / 2)) * ((mplanet) ** (5 / 2)) * ((radiusMoon) ** (5)) *
                   ((eccentricity) ** (2))) / (QualityF * ((semimajor) ** (15 / 2)))
         normed_PowerCoordPlus = []
         normed_PowerCoordPlus.append(semiMajorau)
         normed_PowerCoordPlus.append(distancePlanet)
         normed_PowerCoordPlus.append(np.arcsinh(power+lTidal))
         normed_MoonPowStelFluxPlusTidal.append(power+lTidal)
         normed_MoonPowStelFluxMinusTidal.append(power-lTidal)
         normed_MoonPowStelFluxPlusTidal1D.append(np.arcsinh(power+lTidal)) # #(power+lTidal)/(1.75*10**17)
         #normed_MoonPowStelFluxPlusTidal1D.append((power + lTidal))
         normed_MoonPowStelFlux1D.append(np.arcsinh(power)) #/(1.75*10**17)
         normed_MoonPowStelFluxTidal1D.append(np.arcsinh(lTidal)) #/(1.75*10**17)
         #normed_MoonPowStelFluxPlusTidal1D.append(normed_PowerCoordPlus)
         normed_PowerCoordMinus = []
         lResetPower = (power-lTidal)  #/(1.75*10**17)
         lArcSinhPow =  np.arcsinh(lResetPower) #lResetPower/(1.75*10**17)
         #if (lArcSinhPow < 0):
         #    lArcSinhPow = lArcSinhPow + 60
         normed_PowerCoordMinus.append(lResetPower)
         normed_MoonPowStelFluxMinusTidal1D.append(lArcSinhPow)
         lstTempMoon.append(tempMoon)
         index+=1
     normed_MoonPowStelFluxPlusTidal2D.append(normed_MoonPowStelFluxPlusTidal1D)
     normed_MoonPowStelFluxMinusTidal2D.append(normed_MoonPowStelFluxMinusTidal1D)
     normed_MoonPowStelFlux2D.append(normed_MoonPowStelFlux1D)
     normed_MoonPowTidal2D.append(normed_MoonPowStelFluxTidal1D)

normed_hoststartemp = np.array(lstHostStarPlanetnew)
normed_semimajorMoon = np.array(lstMoon1dnew)
normed_MoonTempStelFlux = np.array(lstTempMoon)
normed_MoonPowStelFluxNP = np.array(normed_MoonPowStelFlux)
normed_MoonPowStelFluxPlusTidalNP = np.array(normed_MoonPowStelFluxPlusTidal)
normed_MoonPowStelFluxMinusTidalNP = np.array(normed_MoonPowStelFluxMinusTidal)
normed_MoonPowStelFluxPlusTidal2DNP = np.array(normed_MoonPowStelFluxPlusTidal2D)
normed_MoonPowStelFluxMinusTidal2DNP = np.array(normed_MoonPowStelFluxMinusTidal2D)
normed_MoonPowStelFlux2DNP = np.array(normed_MoonPowStelFlux2D)
normed_MoonPowTidal2DNP =  np.array(normed_MoonPowTidal2D)

#normed_MoonPowStelFluxPlusTidal2DNP.sort()
#normed_MoonPowStelFluxMinusTidal2DNP.sort()

index = 0
for i in range(178):
    np.append(normed_hoststartemp,[5778])
    np.append(normed_semimajorMoon,[lstSM1D[index]])
    np.append(normed_MoonTempStelFlux,[lstTemp1D[index]])
    index +=1

lsSolarMoonNew = []
index = 0
for i in range(178):
    lsSolarMoonNew.append(lstSM1D[index])
    index +=1

print ("Length of Planet Semi Major", normed_semimajor.shape)
print ("Length of Moon Semi Major", normed_semimajorMoon.shape)
print ("Length of ZLine", normed_MoonPowStelFluxNP.shape)

print ("Power 0 = ",normed_MoonPowStelFluxPlusTidal2DNP[0][0], "Stellar Power = ", normed_MoonPowStelFlux2DNP[0][0], "Tidal Power = ", normed_MoonPowTidal2DNP[0][0], "Distance Planet = ", npSemiMajorPlanet[0], "Distance planet moon", npSemiMajorMoon[0][0])
print ("Power 2 = ",normed_MoonPowStelFluxPlusTidal2DNP[2][0], "Stellar Power = ", normed_MoonPowStelFlux2DNP[2][0], "Tidal Power = ", normed_MoonPowTidal2DNP[2][0], "Distance Planet = ", npSemiMajorPlanet[2], "Distance planet moon", npSemiMajorMoon[2][0])
print ("Power 8 = ",normed_MoonPowStelFluxPlusTidal2DNP[8][0], "Stellar Power = ", normed_MoonPowStelFlux2DNP[8][0], "Tidal Power = ", normed_MoonPowTidal2DNP[8][0], "Distance Planet = ", npSemiMajorPlanet[8], "Distance planet moon", npSemiMajorMoon[8][0])
print ("Power 100 = ",normed_MoonPowStelFluxPlusTidal2DNP[0][99], "Stellar Power = ", normed_MoonPowStelFlux2DNP[0][99], "Tidal Power = ", normed_MoonPowTidal2DNP[0][99], "Distance Planet = ",  npSemiMajorPlanet[0], "Distance planet moon", npSemiMajorMoon[0][99])
print("Power on Earth = ",powerFromAlbedoL(1,6300*1000,0.3))

#Test IO Tidal Heating
K_2OverQ = 0.01
semimajor = 422000*1000
radiusMoon = 1800*1000
mplanet = 1.9*(10**27)
eccentricity = 0.0041
lIOTidal = ((21 / 2) * (K_2OverQ) * ((G) ** (3 / 2)) * ((mplanet) ** (5 / 2)) * ((radiusMoon) ** (5)) *
          ((eccentricity) ** (2))) / (((semimajor) ** (15 / 2)))

print("IO Tidal Heating",lIOTidal)

print ("Power 0 = ",normed_MoonPowStelFluxPlusTidalNP[0], "Distance 0 = ", lstPlanet1dnew[0]+lstMoon1dnew[0])
print ("Power 2 = ",normed_MoonPowStelFluxPlusTidalNP[2], "Distance 1 = ", lstPlanet1dnew[1]+lstMoon1dnew[1])
print ("Power 8 = ",normed_MoonPowStelFluxPlusTidalNP[8], "Distance 2 = ", lstPlanet1dnew[2]+lstMoon1dnew[2])
print ("Power 20 = ",normed_MoonPowStelFluxPlusTidalNP[20], "Distance 20 = ", lstPlanet1dnew[20]+lstMoon1dnew[20])
print("Min Stellar Power = ", min(normed_MoonPowStelFluxNP), "Max stellar power = ", max(normed_MoonPowStelFluxNP))
print("Min Stellar Power + Tidal = ", min(normed_MoonPowStelFluxPlusTidalNP), "Max stellar power plus Tidal = ", max(normed_MoonPowStelFluxPlusTidalNP))

minX = min(normed_semimajorMoon)
maxX = max(normed_semimajorMoon)
minY = min(normed_semimajor)
maxY = max(normed_semimajor)

print("2D Numpy Power Array Plus", normed_MoonPowStelFluxPlusTidal2DNP)
print("2D Numpy Power Array Minus",normed_MoonPowStelFluxMinusTidal2DNP)
print("The shape of Power",normed_MoonPowStelFluxPlusTidal2DNP.shape)

print("Test Power Array Plus",normed_MoonPowStelFluxPlusTidal2DNP)

interp = 'bilinear'
fig, axs = plt.subplots()
plt.axis([0,0.4,0,30])
print("Arc Sin Earth Insolation =",np.arcsinh(1.75*10**17))
points = axs.imshow(normed_MoonPowStelFluxPlusTidal2DNP/np.arcsinh(1.75*10**17), extent = [0,0.4,0,30], origin='lower', interpolation=interp, aspect = "auto")
axs.set_xlabel('Semi Major Axis of Exomoons [AU]')
axs.set_ylabel('Semi Major Axis of Planets [AU]')
print("points = ",points)
cbar = fig.colorbar(points)
cbar.set_label("Earth's insolation by the Sun")
plt.show()
plt.close()


interp = 'bilinear'
fig, axs = plt.subplots()
plt.axis([0,0.4,0,30])
points = axs.imshow(normed_MoonPowStelFluxMinusTidal2DNP/np.arcsinh(lIOTidal), extent = [0,0.4,0,30], origin='lower', interpolation=interp, aspect = "auto")
axs.set_xlabel('Semi Major Axis of Exomoons [AU]')
axs.set_ylabel('Semi Major Axis of Planets [AU]')
cbar = fig.colorbar(points)
cbar.set_label("Tidal Heating received by IO from Jupiter")
plt.show()
print("Test Power Array Minus",normed_MoonPowStelFluxMinusTidal2DNP)

interp = 'bilinear'
fig, axs = plt.subplots()
plt.axis([0,0.4,0,30])
points = axs.imshow(normed_MoonPowStelFluxMinusTidal2DNP/np.arcsinh(lIOTidal), extent = [0,0.4,0,30], origin='lower', interpolation=interp, aspect = "auto")
axs.set_xlabel('Semi Major Axis of Exomoons [AU]')
axs.set_ylabel('Semi Major Axis of Planets [AU]')
cbar = fig.colorbar(points)
cbar.set_label("Tidal Heating received by IO from Jupiter")
plt.show()


interp = 'bilinear'
fig, axs = plt.subplots()
plt.axis([0,0.4,0,30])
points = axs.imshow(normed_MoonPowStelFlux2DNP/np.arcsinh(1.75*10**17), extent = [0,0.4,0,30], origin='lower', interpolation=interp, aspect = "auto")
axs.set_xlabel('Semi Major Axis of Exomoons [AU]')
axs.set_ylabel('Semi Major Axis of Planets [AU]')
print("points = ",points)
cbar = fig.colorbar(points)
cbar.set_label("Earth's insolation by the Sun")
plt.show()
plt.close()


interp = 'bilinear'
fig, axs = plt.subplots()
plt.axis([0,0.4,0,30])
points = axs.imshow(normed_MoonPowTidal2DNP/np.arcsinh(lIOTidal), extent = [0,0.4,0,30], origin='lower', interpolation=interp, aspect = "auto")
axs.set_xlabel('Semi Major Axis of Exomoons [AU]')
axs.set_ylabel('Semi Major Axis of Planets [AU]')
cbar = fig.colorbar(points)
cbar.set_label("Tidal Heating received by IO from Jupiter")
plt.show()

distanceStarMoonAU = random.uniform(1, 1)  # random.uniform(0.2, 30.5)
radiusMoon = random.uniform(1800 * 1000, 1800 * 1000)

lstMBodiies = []
lstEscapeVelBodies = []
lstRootMean = []
colors = []

for i in range(100000):
    bodymass = rndmPLaw(3 * (10 ** 22), 5.9 * (10 ** 25), g=-0.4, size=1)
    mBody = bodymass[0] #random.uniform(3 * (10 ** 19), 5.9 * (10 ** 24)) #, g=-0.5, size=1) #random.uniform(3 * (10 ** 19), 2 * (10 ** 27))
    distanceStarMoon = distanceStarMoonAU*1.496 * ((10) ** (8)) * 1000
    escapeVel = ((2*G*mBody)/(radiusMoon*1000*1000))**(1/2)
    power = powerFromAlbedoL(distanceStarMoon, radiusMoon, 0.3)
    lstMBodiies.append(mBody)
    lstEscapeVelBodies.append(escapeVel)
    print("Mass = ", mBody, "Radius = ",radiusMoon, "Distance From Star", distanceStarMoon, "Escape Vel in km/s:",escapeVel)

    temp = 1000
    #rootMeanVelofNitrogen = ((2*B_Constant*temp)/(28*1.66*10**(-27)))**(1/2)
    #print ("Root Mean Nitrogen m/s = ",rootMeanVelofNitrogen," at", temp)
    #lstRootMean.append(rootMeanVelofNitrogen/1000)
    #colors.append('r')

    rootMeanVelofOxygen = ((2*B_Constant*temp)/(32*1.66*10**(-27)))**(1/2)
    print ("Root Mean Oxygen m/s  ",rootMeanVelofOxygen," at",temp )
    lstRootMean.append(rootMeanVelofOxygen/1000)
    colors.append('g')
    #lstMBodiies.append(mBody)
    #lstEscapeVelBodies.append(escapeVel)

    rootMeanVelofHeleum = ((2*B_Constant*temp)/(4*1.66*10**(-27)))**(1/2)
    print ("Root Mean Helium m/s  = ",rootMeanVelofHeleum," at",temp )
    lstRootMean.append(rootMeanVelofHeleum/1000)
    colors.append('red')
    lstMBodiies.append(mBody)
    lstEscapeVelBodies.append(escapeVel)

    rootMeanVelofHydrogen = ((2*B_Constant*temp)/(2*1.66*10**(-27)))**(1/2)
    print ("Root Mean Hydrogen m/s = ",rootMeanVelofHydrogen," at",temp)
    lstRootMean.append(rootMeanVelofHydrogen/1000)
    colors.append('y')
    lstMBodiies.append(mBody)
    lstEscapeVelBodies.append(escapeVel)

    #rootMeanVelofCarbonDiOxide = ((2*B_Constant*temp)/(44*1.66*10**(-27)))**(1/2)
    #print ("Root Mean Carbon Dioxide m/s = ",rootMeanVelofCarbonDiOxide," at",temp)
    #lstRootMean.append(rootMeanVelofCarbonDiOxide/1000)
    #lstTemp.append(temp)
    #colors.append('black')
    #lstMBodiies.append(mBody)
    #lstEscapeVelBodies.append(escapeVel)

npMBodies = np.array(lstMBodiies)
npEscapeVel = np.array(lstEscapeVelBodies)
npRootMean = np.array(lstRootMean)

fig,ax = plt.subplots()

plt.scatter(npMBodies/(5.9*10**(24)), npEscapeVel, s=1, c='b')
plt.scatter(npMBodies/(5.9*10**(24)), npRootMean*6, s=1, c=colors)

red_patch = mpatches.Patch(color = 'red', label = 'Mean Velocity of Helium')
yellow_patch = mpatches.Patch(color = 'yellow', label = 'Mean Velocity of Hydrogen')
green_patch = mpatches.Patch(color = 'green', label = 'Mean Velocity of Oxygen')
blue_patch = mpatches.Patch(color = 'blue', label = 'Escape Velocity of body')
ax.legend(handles = [red_patch,yellow_patch,green_patch,blue_patch])

#plt.plot(npMBodies/(5.9*10**(24)), npEscapeVel, c='b')
xs = npMBodies  # semiMajor
ys = npEscapeVel  # temperatureMoons

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

plt.xlabel("Mass of the bodies [M_E]")
plt.ylabel("Escape Velocity [Km/S]")
plt.show()

#mEarth = random.uniform(5.9 * (10 ** 24), 5.9 * (10 ** 24))
#mJupiter = random.uniform(2 * (10 ** 27), 2 * (10 ** 27))
mMoon = random.uniform(7.34 * (10 ** 22), 7.34 * (10 ** 22))
#mMars = random.uniform(6.39 * (10 ** 23), 6.39 * (10 ** 23))
mEarth = mMoon

#distanceStarEarth = 1
#distanceStarJup = 5.2
distanceStarMoon = 1
#distanceStarMars = 1.5

distanceStarEarth = distanceStarMoon
radiusEarth = 6370*1000
radiusJup = 71000*1000
radiusMoon = 1800*1000
radiusMars = 3830*1000
radiusEarth = radiusMoon

escapeVel = ((2 * G * mEarth) / (radiusEarth * 1000 * 1000)) ** (1 / 2)
power = powerFromAlbedoL(distanceStarEarth, radiusEarth, 0.3)
print("Mass = ", mEarth, "Radius = ", radiusEarth, "Distance From Star", distanceStarEarth, "Escape Vel in km/s:", escapeVel, "Stellar Power=",power)
escapeVel = escapeVel * 1000  # Escape Vel in M/s

lstRootMean = []
lstTemp = []
lstEscapeRate = []
colors = []

numberDensity = 10**8 #Number of Molecules in exosphere/m^3

for temp in range (200,2200,1):
    #Root Mean Most Probable Vel of Gases in Exosphere @ 300k
    #rootMeanVelofNitrogen = ((2*B_Constant*temp)/(28*1.66*10**(-27)))**(1/2)
    #print ("Root Mean Nitrogen m/s = ",rootMeanVelofNitrogen," at", temp)
    #lstRootMean.append(rootMeanVelofNitrogen)
    #lstTemp.append(temp)
    #colors.append('r')

    rootMeanVelofOxygen = ((2*B_Constant*temp)/(32*1.66*10**(-27)))**(1/2)
    print ("Root Mean Oxygen m/s  ",rootMeanVelofOxygen," at",temp )
    lstRootMean.append(rootMeanVelofOxygen)
    lstTemp.append(temp)
    colors.append('g')

    rootMeanVelofHelium = ((2*B_Constant*temp)/(4*1.66*10**(-27)))**(1/2)
    print ("Root Mean Helium m/s  = ",rootMeanVelofHelium," at",temp )
    lstRootMean.append(rootMeanVelofHelium)
    lstTemp.append(temp)
    colors.append('r')

    rootMeanVelofHydrogen = ((2*B_Constant*temp)/(2*1.66*10**(-27)))**(1/2)
    print ("Root Mean Hydrogen m/s = ",rootMeanVelofHydrogen," at",temp)
    lstRootMean.append(rootMeanVelofHydrogen)
    colors.append('y')
    lstTemp.append(temp)

    #rootMeanVelofCarbonDiOxide = ((2*B_Constant*temp)/(44*1.66*10**(-27)))**(1/2)
    #print ("Root Mean Carbon Dioxide m/s = ",rootMeanVelofCarbonDiOxide," at",temp)
    #lstRootMean.append(rootMeanVelofCarbonDiOxide)
    #lstTemp.append(temp)
    #colors.append('black')


    PhiEscape = (numberDensity * rootMeanVelofOxygen) * ((escapeVel ** 2 / rootMeanVelofOxygen ** 2) + 1) * np.exp(
        -escapeVel ** 2 / rootMeanVelofOxygen ** 2)
    print ("Oxygen Escape Rate is /m^2/s", PhiEscape, "at", temp, "Escape Vel of Body = ",escapeVel)
    lstEscapeRate.append(PhiEscape)

    #PhiEscape = (numberDensity * rootMeanVelofNitrogen) * ((escapeVel ** 2 / rootMeanVelofNitrogen ** 2) + 1) * np.exp(
    #    -escapeVel ** 2 / rootMeanVelofNitrogen ** 2)
    #print ("Nitrogen Escape Rate is /m^2/s", PhiEscape, "at", temp)
    #lstEscapeRate.append(PhiEscape)

    PhiEscape = (numberDensity * rootMeanVelofHelium) * ((escapeVel ** 2 / rootMeanVelofHelium ** 2) + 1) * np.exp(
        -escapeVel ** 2 / rootMeanVelofHelium ** 2)
    print ("Helium Escape Rate is /m^2/s", PhiEscape, "at", temp)
    lstEscapeRate.append(PhiEscape)

    PhiEscape = (numberDensity * rootMeanVelofHydrogen) * ((escapeVel ** 2 / rootMeanVelofHydrogen ** 2) + 1) * np.exp(
        -escapeVel ** 2 / rootMeanVelofHydrogen ** 2)
    print ("Hydrogen Escape Rate is /m^2/s", PhiEscape, "at", temp)
    lstEscapeRate.append(PhiEscape)

    #PhiEscape = (numberDensity * rootMeanVelofCarbonDiOxide) * ((escapeVel ** 2 / rootMeanVelofCarbonDiOxide ** 2) + 1) * np.exp(
    #    -(escapeVel) ** 2 / rootMeanVelofCarbonDiOxide ** 2)
    #print ("Carbondioxide Escape Rate is /m^2/s", PhiEscape, "at", temp)
    #lstEscapeRate.append(PhiEscape)

npRootMean = np.array(lstRootMean)
npTemp = np.array(lstTemp)
npEscapeRate = np.array(lstEscapeRate)

fig,ax = plt.subplots()
plt.scatter(npTemp, npRootMean, s=1, c=colors)

red_patch = mpatches.Patch(color = 'red', label = 'Mean Velocity of Helium')
yellow_patch = mpatches.Patch(color = 'yellow', label = 'Mean Velocity of Hydrogen')
green_patch = mpatches.Patch(color = 'green', label = 'Mean Velocity of Oxygen')
ax.legend(handles = [red_patch,yellow_patch,green_patch])

xs = npMBodies  # semiMajor
ys = npEscapeVel  # temperatureMoons

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

plt.xlabel("Temperature of Exospheres [K]")
plt.ylabel("Molecular Velocity [M/S]")
plt.show()


fig,ax =plt.subplots()
plt.scatter(npTemp, np.log(npEscapeRate),s=1, c=colors)

xs = npTemp  # semiMajor
ys = npEscapeRate  # temperatureMoons

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()

plt.xlim(xmin, xmax)
plt.ylim(ymin, ymax)

red_patch = mpatches.Patch(color = 'red', label = 'Escape Rate of Helium')
yellow_patch = mpatches.Patch(color = 'yellow', label = 'Escape Rate of Hydrogen')
green_patch = mpatches.Patch(color = 'green', label = 'Escape Rate of Oxygen')
ax.legend(handles = [red_patch,yellow_patch,green_patch])

plt.xlabel("Temperature of Exospheres [K]")
plt.ylabel("Escape Rate in log scale[Particles /M^2/S]")
plt.show()