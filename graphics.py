import matplotlib.pyplot as plt
import matplotlib.ticker 
import numpy as np
from calculations import *
import os

if not os.path.isdir("grafic_img"): # папка для хранения графиков
    os.mkdir("grafic_img")

velocitiyScalarArray = [lengthOfVector(velocityArray[i]) for i in range(period + 1)]
# Данные по спискам из программы расчётов
time = np.array(timeArray)
height = np.array(heightArray)
mass = np.array(massOfRocketArray)
velocity = np.array(velocitiyScalarArray)
acceleration = np.array(accelerationArray)
angle = np.array(angelArray)
pressure = np.array(pressureArray)
ro = np.array(roArray)
dragForce = np.array(dragForceArray)
position = np.array(positionArray)
forceOfGravity = np.array(forceOfGravityArray)
ForceOFThrust = np.array(forceOfThrustArray)

# Построение графиков функции
height_from_time = plt.figure("Высота", figsize=(10, 5), dpi=150) # номер фиуры

plt.plot(time, height, 'g') # настройки отображения: цвет, линия

# названия осей и графика
plt.title("Высота от времени")
plt.xlabel('Время, с')
plt.ylabel('Высота, м')

height_from_time.savefig("grafic_img/height_from_time.png")

mass_from_time = plt.figure("Масса", figsize=(10, 5), dpi=150)

plt.plot(time, mass, 'r') # настройки отображения: цвет, линия

# названия осей и графика
plt.title("Масса от времени")
plt.xlabel('Время, с')
plt.ylabel('Масса, кг')
mass_from_time.savefig("grafic_img/mass_from_time.png")

pressure_from_height = plt.figure("Давление", figsize=(10, 5), dpi=150)

plt.plot(height, pressure, 'blue') # настройки отображения: цвет, линия

# названия осей
plt.title("Давление от высоты")
plt.xlabel('Высота, м')
plt.ylabel('Давление, Па')
pressure_from_height.savefig("grafic_img/pressure_from_height.png")

dragforse_from_height = plt.figure("Сила сопротивления", figsize=(10, 5), dpi=150)

plt.plot(height, dragForce, 'orange') # настройки отображения: цвет, линия

# названия осей
plt.title("Сила сопросивления от высоты")
plt.xlabel('Высота, м')
plt.ylabel('Сила сопротивления, Н')
dragforse_from_height.savefig("grafic_img/dragforse_from_height")

velocity_from_time = plt.figure('velocity', figsize=(10, 5), dpi=150)

plt.plot(time, velocity, 'purple') # настройки отображения: цвет, линия

# названия осей
plt.title("Скорость от времени")
plt.xlabel('Время, с')
plt.ylabel('Скорость, м/с')
velocity_from_time.savefig("grafic_img/velocity.png")

forceOfGravity_from_time = plt.figure("Сила тяготения", figsize=(8, 5), dpi=150)

plt.plot(time, forceOfGravity, 'blue') # настройки отображения: цвет, линия

# названия осей
plt.title("Сила тяготения от времени")
plt.xlabel('Время, с')
plt.ylabel('Сила тяготения, Н')
forceOfGravity_from_time.savefig("grafic_img/forceOfGravity_from_time")

ro_from_height = plt.figure("Плотность", figsize=(10, 5), dpi=150)

plt.plot(height, ro, 'm') # настройки отображения: цвет, линия

# названия осей
plt.title("Плотность воздуха от высоты")
plt.xlabel('Высота, м')
plt.ylabel('Плотность, кг/м^3')
ro_from_height.savefig("grafic_img/ro_from_height")

force_of_thrust = plt.figure("Силя тяги", figsize=(10, 5), dpi=150)

plt.plot(time, ForceOFThrust, 'olive') # настройки отображения: цвет, линия

# названия осей
plt.title("Сила тяги двигателя от времени")
plt.xlabel('Время, с')
plt.ylabel('Сила тяги двигателя, Н')
force_of_thrust.savefig("grafic_img/force_of_thrust")

plt.show()