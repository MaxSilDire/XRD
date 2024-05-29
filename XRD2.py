# -*- coding: utf-8 -*-

import numpy
import xrayutilities as xru
from xrayutilities.materials.cif import CIFFile
from xrayutilities.materials.material import Crystal
from statistics import mean
import matplotlib.pyplot as plt 
import matplotlib 
from multiprocessing import freeze_support
plt.style.use('ggplot')
matplotlib.rcParams.update({'font.size': 18})
fig_size = [15, 12]
plt.rcParams["figure.figsize"] = fig_size

def main():
    #симуляция для кубической структуры
    xu_cif = CIFFile("Si Cubic MP-149.cif")
    xu_crystal = Crystal(name="Si Cubic", lat=xu_cif.SGLattice())
    two_theta = numpy.arange(30, 140.04, 0.02)
    powder = xru.simpack.smaterials.Powder(xu_crystal, 1)
    pm = xru.simpack.PowderModel(powder, I0=35000)    
    intensitiesSim1 = pm.simulate(two_theta)       
    pm.close()
    
    # отсеиваем слабые пики
    Sim1 = [x for x in intensitiesSim1 if x <= 10]
    
    # считаем площадь спектра для первой структуры, интегрируя вдоль оси
    AreaSim1 = numpy.trapz(Sim1, dx=0.02)
    
    #симуляция для второй структуры
    xu_cif1 = CIFFile("Si Tetragonal MP-92.cif")
    xu_crystal1 = Crystal(name="Si", lat=xu_cif1.SGLattice())
    two_theta = numpy.arange(30, 140.04, 0.02)
    powder = xru.simpack.smaterials.Powder(xu_crystal1, 1)
    pm = xru.simpack.PowderModel(powder, I0=35000)    
    intensitiesSim2 = pm.simulate(two_theta)       
    pm.close()
    
    # построение графика двух симуляций для наглядного сравнения
    xru.simpack.plot_powder(two_theta, intensitiesSim1, intensitiesSim2, scale='sqrt', show_diff=False)
   
    # отсеиваем слабые пики
    Sim2 = [x for x in intensitiesSim2 if x <= 10]
    
    # считаем площадь спектра для второй структуры, интегрируя вдоль оси
    AreaSim2 = numpy.trapz(Sim2, dx=0.02)
    
    # определяем спектр с большей площадью для вычисления геометрической вероятности
    if AreaSim1 > AreaSim2:
        AreaMax = AreaSim1
        AreaMin = AreaSim2
    elif AreaSim1 < AreaSim2:
        AreaMax = AreaSim2
        AreaMin = AreaSim1
    
    # вычисляем процент возможной ошибки
    P = (AreaMin / AreaMax) * 100
    print(AreaMin)
    print(AreaMax)
    print(P)
if __name__ == '__main__':
    freeze_support()
    main()

    
