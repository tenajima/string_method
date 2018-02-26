import numpy as np
import time
import copy
import matplotlib.pyplot as plt
import os
import shutil
from simtk.unit import *
from simtk.openmm.app import *
from simtk.openmm import *
from scipy.interpolate import interp1d


class String:
    def __init__(self, Points, restart=False, count=0, work_dir=''):
        self.restart = restart
        self.work_dir = work_dir
        self.Points = []
        self.centerPoints = []
        self.log_of_convergence = []
        self.diff_end_point = []
        self.count = count
        self.Points = []
        self.centerPoints = []
        self.log_of_convergence = []
        self.diff_end_point = []
        self.numberOfPoints = len(Points)
        print(self.numberOfPoints)
        self.check_log_dir()
        self.solute_sum = len(Points[0].solute_array)
        # self.dt = 0.05 * min(0.2, self.numberOfPoints ** (-1))
        self.dt = 0.00001
        print('dt is {}'.format(self.dt))
        self.Points = Points
        # self.logOfCenterPoints()

    def stringMethod(self):
        d = 10000
        min_d = 100000
        oscillation_count = 0
        tol = max(10 ** (-10), self.numberOfPoints ** (-4.0))
        print("tol is " + str(tol))
        com_log = open('com_log.dat', 'w')
        while True:
            self.count += 1
            difference = []
            # oldPoints = copy.deepcopy(self.Points)
            old_solute_array = [self.Points[i].get_solute_array() for i in range(number_of_points)]
            for i in range(self.numberOfPoints):
                self.Points[i].time_development(self.dt, in_water=False, rism_prefix='ADP_string')
                # self.logOfCenterPoints()
                # self.plotString(count)
            # self.plotString()
            self.interpolate()

            com_x, com_y, com_z = self.logOfCenterPoints()
            com_log.write(f'{self.count:>4d}  {com_x:>15.3e}  {com_y:>15.3e}  {com_z:>15.3e}\n')


            for i in range(self.numberOfPoints):
                diff = 0.0
                for solute_num in range(self.solute_sum):
                    for coordinate in range(3):
                        diff += (self.Points[i].solute_array[solute_num][coordinate]._value - old_solute_array[i][solute_num][coordinate]) ** 2.0
                difference.append(diff ** 0.50)
            d = max(difference) / self.dt
            self.diff_end_point.append(difference[0] + difference[self.numberOfPoints - 1])
            self.log_of_convergence.append(d)
            print(f'count is {self.count} : d is {d:10.3e}')
            if d > min_d:
                oscillation_count += 1
                print('oscillation count is {}'.format(oscillation_count))
                if oscillation_count > 10:
                    print("d > old d! and d is " + str(d))
                    print("count is " + str(self.count))
                    self.makeLogFile(count=self.count)
                    # self.makeLogLastFile(count)
                    # self.plotCenterPoints()
                    break
            elif d < tol:
                print("d < TOL. And d is " + str(d))
                print("count is " + str(self.count))
                # self.makeLogLastFile(count)
                # self.interpolate()
                self.makeLogFile(count=self.count)
                # self.plotCenterPoints()
                break
            else:
                # oscillation_count=0
                if self.count % 100 == 0:
                    self.makeLogFile(count=self.count)
                if d < min_d:
                    min_d = d
                    oscillation_count=0

            # if count > 1000:
            #     self.makeLogFile(count=count)
            #     break

    def interpolate(self):
        s = [0.0]
        alpha = []
        for i in range(1, self.numberOfPoints):
            length = 0.0
            for solute_num in range(self.solute_sum):
                for coordinate in range(3):
                    length += (self.Points[i].solute_array[solute_num][coordinate]._value -
                               self.Points[i - 1].solute_array[solute_num][coordinate]._value) ** 2.0
            s.append(length ** 0.50 + s[i - 1])
        for i in range(self.numberOfPoints):
            alpha.append(s[i] / s[self.numberOfPoints - 1])

        for solute_num in range(self.solute_sum):
            for coordinate in range(3):
                debug_plot = []
                debug_x1 = []
                debug_plot2 = []
                debug_x2 = []
                measure = []
                for i in range(self.numberOfPoints):
                    measure.append(self.Points[i].solute_array[solute_num][coordinate]._value)
                # cubicInterp = interp1d(alpha, measure)
                # cubicInterp = interp1d(alpha, measure, kind='quadratic')
                cubicInterp = interp1d(alpha, measure, kind='cubic')

                for i in range(self.numberOfPoints):
                    self.Points[i].solute_array[solute_num][coordinate] = Quantity(value=cubicInterp(i / (self.numberOfPoints - 1)),
                                                                                   unit=nanometer)
                    # debug_x1.append(i/(self.numberOfPoints - 1))
                    # debug_plot.append(cubicInterp(i/(self.numberOfPoints - 1)))
                    # debug_x2.append(alpha[i])
                    # debug_plot2.append(cubicInterp(alpha[i]))
                    # print(f'alpha[i] is                    {alpha[i]}')
                    # print(f'i/(self.numberOfPoints - 1) is {i/(self.numberOfPoints - 1)}')
                # if coordinate == 0 and solute_num == 0:
                # plt.plot(debug_x1, debug_plot, '-o', color='red', label='i/(self.numberOfPoints - 1)')
                # plt.plot(debug_x2, debug_plot2, '-o', color='blue', label='alpha[i]')
                # plt.legend()
                # plt.show()
        for i in range(self.numberOfPoints):
            self.Points[i].update_position()

#------------------------------------------- debug -------------------------------------------------------
        # s = [0.0]
        # for i in range(1, self.numberOfPoints):
        #     length = 0.0
        #     for solute_num in range(self.solute_sum):
        #         for coordinate in range(3):
        #             length += (self.Points[i].solute_array[solute_num][coordinate]._value -
        #                        self.Points[i - 1].solute_array[solute_num][coordinate]._value) ** 2.0
        #     s.append(length ** 0.50)
        # print('after interpolate: ', s)
        # input()
#---------------------------------------------------------------------------------------------------------


    def plotString(self, count=0):
        import matplotlib.pyplot as plt
        x = []
        y = []
        for i in range(self.numberOfPoints):
            length = np.zeros(3)
            for coordinate in range(3):
                length[0] += (self.Points[i].solute_array[0][coordinate] - self.Points[i].solute_array[1][coordinate]) ** 2.0
                length[1] += (self.Points[i].solute_array[1][coordinate] - self.Points[i].solute_array[2][coordinate]) ** 2.0
                length[2] += (self.Points[i].solute_array[0][coordinate] - self.Points[i].solute_array[2][coordinate]) ** 2.0
            length **= 0.5
            x.append(length[0])
            y.append(length[1])
        plt.xlim(1.5, 5.0)
        plt.ylim(1.5, 5.0)
        plt.gca().set_aspect('equal', adjustable='box')
        plt.plot(x, y, "o")
        plt.plot(x, y)
        # plt.savefig("img/string_%03d.png"%count)
        plt.show()
        plt.clf()

    def check_log_dir(self):
        debug = True
        self.top_dir_name = f'{self.work_dir}point_num_{self.numberOfPoints}/'
        if self.restart:
            return
        if os.path.exists(self.top_dir_name):
            print('そのストリングのデータはあります.削除しますか?(新しく作り直すなら:yes)')
            if debug == False:
                user_input = input()
            else:
                user_input = 'yes'

            if user_input == 'yes':
                self.make_new_log = True
                shutil.rmtree(self.top_dir_name)
                os.mkdir(self.top_dir_name)
            else:
                self.make_new_log = False
                return
        else:
            if debug:
                self.make_new_log = True
            print(os.getcwd())
            os.mkdir(self.top_dir_name)

    def makeLogFile(self, count=0):
        if self.restart:
            for i in range(self.numberOfPoints):
                filename = (f'solute_num_{i:02d}.pdb')
                self.Points[i].export_pdb(dir_name=self.top_dir_name, file_name=filename)
                return
        if self.make_new_log == False:
            print('logファイルは作り直しませんでした')
            return
        file_array = []
        # top_dir_name = f'./point_num_{self.numberOfPoints}/'
        dir_name = f'{self.top_dir_name}data_{count}/'
        os.mkdir(dir_name)
        # os.chdir(dir_name)

        config_file = open(dir_name + 'config.txt', 'w')
        config_file.write('calculation config\n')
        config_file.write(f'Number of point: {self.numberOfPoints}\n')

        for i in range(self.numberOfPoints):
            # filename = ("solute_num" + str(count) + ".3d")
            filename = (f'solute_num_{i:03d}.pdb')
            self.Points[i].export_pdb(dir_name=dir_name, file_name=filename)
        # os.chdir('../../')


        # f = open(filename, "w")
        # for i in range(self.numberOfPoints):
        #     f.write("{0}  {1}  {2}".format(self.Points[i].r_Cl1[0], self.Points[i].r_CH3[0], self.Points[i].r_Cl2[0]))
        #     f.write("\n")
        # f.close()

    def makeLogLastFile(self, count):
        filename = "logLastPhi" + str(count) + ".dat"
        f = open(filename, "w")
        for i in range(self.numberOfPoints):
            f.write("{0}  {1}  {2}".format(self.Points[i].r_Cl1[0], self.Points[i].r_CH3[0], self.Points[i].r_Cl2[0]))
            f.write("\n")
        f.close()
    def logOfCenterPoints(self):
        average_of_center_of_mass = []
        for point in self.Points:
            average_of_center_of_mass.append(point.get_center_of_mass())
        average_of_center_of_mass = np.array(average_of_center_of_mass)

        return (np.average(average_of_center_of_mass[:,0]), np.average(average_of_center_of_mass[:,1]),
                np.average(average_of_center_of_mass[:,2]))
    def plotCenterPoints(self):
        import matplotlib.pyplot as plt
        plt.plot(self.centerPoints)
        # plt.ylim(-1e-13,1e-13)
        # plt.ylim(-0.05, 0.05)
        # plt.savefig("centerOfreplica.png")
        plt.show()

def log_number_of_string_point_and():
    f = open("variour_with_number_of_points.dat", "w")
    for i in range(15, 30):
        if i % 2 == 1:
            start = time.time()
            string = String(i)
            print("now calculate , the number of string points is {}".format(i))
            string.stringMethod()
            elapsed_time = time.time() - start
            f.write(" {0}  {1}  {2}\n".format(i, string.count, elapsed_time))
        # string.plotString(1)
    f.close()

def plot_step_and_convergence():
    for i in range(15, 20):
        if i % 2 == 1:
            string = String(i)
            print('now calculate, the number of string point is {}'.format(i))
            string.stringMethod()
            plt.plot(string.log_of_convergence, label='number of point {}'.format(string.numberOfPoints))
    plt.xlim(0,100)
    plt.ylim(0,50)
    plt.legend()
    plt.savefig('step_and_convergence.png')
    plt.show()

def plot_end_convergence():
    string = String(15)
    string.stringMethod()
    plt.plot(string.diff_end_point, label = 'number of point {}'.format(string.numberOfPoints))
    plt.legend()
    plt.savefig('convergence_of_end_point.png')
    plt.show()

if __name__ == '__main__':
    print('start string method')
    from ADP import ADP
    from post_slack import Slack
    start_time = time.time()
    points = []
    restart = False
    # os.chdir('./C7eq_C7ax_v10')
    # os.chdir(dir_name)
    number_of_points = 256
    dir_name = f'./vacuum_{number_of_points}/'


    print('pdbファイルを読み込み中')
    if restart:
        count = 600
        # os.chdir(f'./point_num_128/data_{count}')
        dir_name = f'./vacuum_{number_of_points}/point_num_{number_of_points}/data_{count}/'
        for i in range(number_of_points):
            print(f'restart: {i}')
            point = ADP(dir_name=dir_name, number_of_point=i, restart=restart)
            points.append(point)
        dir_name= f'./vacuum_{number_of_points}/'
    else:
        count = 0
        for i in range(number_of_points):
            print(f'initial start: {i}')
            point = ADP(number_of_point=i, dir_name=dir_name)
            points.append(point)
    print('count is :', count)
    string = String(points, restart=restart, count=count, work_dir=dir_name)
    f=open('toSlack.txt','w')
    f.write('string法の計算が終わりました')
    f.close()
    # test = string.Points[0].get_solute_array()
    # print(test)
    # print(string.Points[0].solute_array)
    string.stringMethod()
    # slack = Slack()
    # slack.read_file()
    # slack.post()

    elapsed_time = time.time() - start_time
    print('elapsed time is {} [s]'.format(elapsed_time))