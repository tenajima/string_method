import numpy as np
import time
import os
from abc import ABCMeta, abstractmethod


class string(metaclass=ABCMeta):
    def __init__(self, Points, restart=False, count=0, work_dir=''):
        self.work_dir = work_dir
        self.Points = []
        self.centerPoints = []
        self.log_of_convergence = []
        self.count = count
        self.Points = []
        self.centerPoints = []
        self.log_of_convergence = []
        self.numberOfPoints = len(Points)
        self.check_log_dir()
        self.solute_sum = len(Points[0].solute_array)
        # self.dt = 0.05 * min(0.2, self.numberOfPoints ** (-1))
        self.dt = 0.00001
        print('dt is {}'.format(self.dt))
        self.Points = Points

    def stringMethod(self):
        min_d = 100000
        oscillation_count = 0
        tol = max(10 ** (-10), self.numberOfPoints ** (-4.0))
        print("tol is " + str(tol))
        while True:
            self.count += 1
            difference = []
            # old_solute_array = [self.Points[i].get_solute_array() for i in range(number_of_points)]
            old_solute_array = np.array([self.Points[i].solute_array for i in range(number_of_points)])
            for i in range(self.numberOfPoints):
                self.Points[i].time_development(self.dt, in_water=False, rism_prefix='ADP_string')
            self.interpolate()

            for i in range(self.numberOfPoints):
                diff = 0.0
                for solute_num in range(self.solute_sum):
                    for coordinate in range(3):
                        diff += (self.Points[i].solute_array[solute_num][coordinate]._value - old_solute_array[i][solute_num][coordinate]) ** 2.0
                difference.append(diff ** 0.50)
            d = max(difference) / self.dt
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

    @abstractmethod
    def interpolate(self):
        raise NotImplementedError


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
            filename = (f'solute_num_{i:03d}.pdb')
            self.Points[i].export_pdb(dir_name=dir_name, file_name=filename)


    def makeLogLastFile(self, count):
        filename = "logLastPhi" + str(count) + ".dat"
        f = open(filename, "w")
        for i in range(self.numberOfPoints):
            f.write("{0}  {1}  {2}".format(self.Points[i].r_Cl1[0], self.Points[i].r_CH3[0], self.Points[i].r_Cl2[0]))
            f.write("\n")
        f.close()



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