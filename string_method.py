import numpy as np
import time
import os
import shutil
from scipy.interpolate import interp1d

class String(object):
    def __init__(self, points, work_dir='', count=0):
        self.work_dir = work_dir
        self.points = []
        self.Points = []
        self.centerPoints = []
        self.count = count
        self.number_of_points = len(points)
        self.check_log_dir()
        self.solute_sum = len(points[0].solute_array)
        self.dt = 0.05 * min(0.2, self.number_of_points ** (-1))
        # self.dt = 0.01
        print('dt is {}'.format(self.dt))
        self.Points = points

    def stringMethod(self):
        min_d = 100000
        oscillation_count = 0
        tol = max(10 ** (-10), self.number_of_points ** (-4.0))
        print("tol is " + str(tol))
        while True:
            self.count += 1
            difference = []
            old_solute_arrays = np.array([self.Points[i].solute_array for i in range(self.number_of_points)])
            for i in range(self.number_of_points):
                self.Points[i].time_development(self.dt)
            self.interpolate()

            for i in range(self.number_of_points):
                diff = 0.0
                for solute_num in range(self.solute_sum):
                    for coordinate in range(3):
                        diff += (self.Points[i].solute_array[solute_num][coordinate] - old_solute_arrays[i][solute_num][coordinate]) ** 2.0
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
                    break
            elif d < tol:
                print("d < TOL. And d is " + str(d))
                print("count is " + str(self.count))
                self.makeLogFile(count=self.count)
                break
            else:
                if self.count % 100 == 0:
                    self.makeLogFile(count=self.count)
                if d < min_d:
                    min_d = d
                    oscillation_count=0

    def interpolate(self):
        s = [0.0]
        alpha = []
        for i in range(1, self.number_of_points):
            length = 0.0
            for solute_num in range(self.solute_sum):
                for coordinate in range(3):
                    length += (self.Points[i].solute_array[solute_num][coordinate] -
                               self.Points[i - 1].solute_array[solute_num][coordinate]) ** 2.0
            s.append(length ** 0.50 + s[i - 1])
        for i in range(self.number_of_points):
            alpha.append(s[i] / s[self.number_of_points - 1])

        for solute_num in range(self.solute_sum):
            for coordinate in range(3):
                measure = []
                for i in range(self.number_of_points):
                    measure.append(self.Points[i].solute_array[solute_num][coordinate])
                cubicInterp = interp1d(alpha, measure, kind='cubic')

                for i in range(self.number_of_points):
                    self.Points[i].solute_array[solute_num][coordinate] = cubicInterp(i / (self.number_of_points - 1))


    def check_log_dir(self):
        debug = False
        self.make_new_log = True
        self.top_dir_name = f'{self.work_dir}point_num_{self.number_of_points}/'
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
        if self.make_new_log == False:
            print('logファイルは作り直しませんでした')
            return
        dir_name = f'{self.top_dir_name}data_{count}/'
        os.mkdir(dir_name)
        config_file = open(dir_name + 'config.txt', 'w')
        config_file.write('calculation config\n')
        config_file.write(f'Number of point: {self.number_of_points}\n')
        for i in range(self.number_of_points):
            filename = (f'solute_num_{i:02d}.pdb')
            self.Points[i].export_solute_array(file_name=f'{self.top_dir_name}data_{count}/{filename}')





if __name__ == '__main__':
    start_time = time.time()

    elapsed_time = time.time() - start_time
    print('elapsed time is {} [s]'.format(elapsed_time))