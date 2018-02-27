from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
import shutil
import subprocess
import time

class ADP:
    '''
    アラニンジペプチドのモデル
    '''
    def __init__(self, dir_name='', file_name=None, number_of_point=0, restart=False):
        '''

        :param number_of_point: ストリングの何番目に位置するか
        :param restart: リスタートファイルを使うか
        '''
        # os.chdir('./C7eq_C7ax_initial_string')
        if file_name is not None:
            self.file_name = f'{dir_name}{file_name}'
        elif restart:
            self.file_name = f'{dir_name}solute_num_{number_of_point:03d}.pdb'
        else:
            self.file_name = f'{dir_name}point_{number_of_point:03d}_opt_com.pdb'
            # self.file_name = f'point_{number_of_point}_optimized.pdb'
        self.pdb = PDBFile(self.file_name)
        self.forcefield = ForceField('amber99sb.xml')
        # self.forcefield = ForceField('amber99sbildn.xml')
        self.system = self.forcefield.createSystem(self.pdb.topology, removeCMMotion=True)
        self.integrator = VerletIntegrator(0.001)
        self.simulation = Simulation(self.pdb.topology, self.system, self.integrator)
        self.simulation.context.setPositions(self.pdb.positions)
        self._solute_array = self.pdb.getPositions(asNumpy=True)
        self.mass_list = [12.01078, 15.99943, 12.0107, 1.007947, 1.007947, 1.007947, 14.00672, 12.01078, 12.01078,
                          15.99943, 12.01078, 1.007947, 1.007947, 1.007947, 1.007947, 1.007947,14.00672,
                          12.01078, 1.007947, 1.007947, 1.007947, 1.007947]
        self.number_of_particles = 22 #アラニンジペプチドの粒子数
        self.non_bonded_force = self.system.getForce(3)
        self.inverse_sum_of_mass = 1 / np.sum(self.mass_list)
        self.count = 0
    def run_onestep(self):
        print(self.simulation.context.getState(getForces=True).getForces(asNumpy=True))
        self.simulation.reporters.append(PDBReporter(f'{self.file_name[:-4]}_one_step.pdb', 1))
        print('---------------------')
        print(self.simulation.context.getState(getForces=True).getForces(asNumpy=True))
        self.simulation.step(1)


    def get_potential_energy(self):
        return self.simulation.context.getState(getEnergy=True).getPotentialEnergy().in_units_of(kilocalorie_per_mole)._value

    # def get_solute_array(self,unit_=nanometer):
    #     '''
    #     getter of solute_array
    #     solute_arrayの単位を除いた値を返す.単位系はnmで返すようにする.
    #     :return:
    #     '''
    #     return self.solute_array.in_units_of(unit_)._value
    @property
    def solute_array(self, with_unit=False):
        if with_unit:
            return self._solute_array
        else:
            return self._solute_array._value

    @solute_array.setter
    def solute_array(self,value):
        self._solute_array = Quantity(value=value, unit=nanometer)


    def update_position(self):
        '''
        self.solute_arrayの変更をOpenMMのpositionにも教える.
        :return:
        '''
        position = []
        for i in range(len(self._solute_array)):
            position.append(tuple(self._solute_array[i].in_units_of(nanometer)))
        self.simulation.context.setPositions(position)

    def time_development(self, dt=0.0005):
        '''
        string法において点を1step時間発展させる
        :param dt:
        :return:
        '''
        force_vacuum = self.simulation.context.getState(getForces=True).getForces(asNumpy=True)
        force_vacuum = force_vacuum.in_units_of(kilocalorie_per_mole/angstrom)
        self._solute_array += Quantity(force_vacuum._value * dt, unit=angstrom)
        self.update_position()
        self._solute_array -= Quantity(self.get_center_of_mass(), unit=angstrom)
        self.count += 1

    def get_center_of_mass(self):
        '''
        アラニンジペプチドの重心を計算
        :return:
        '''
        center_of_mass = np.zeros(3)
        for i in range(22):
            for coordinate in range(3):
                center_of_mass[coordinate] += self._solute_array._value[i][coordinate] * self.mass_list[i]
        return center_of_mass * self.inverse_sum_of_mass

    def rearange_center_of_mass(self):
        com = self.get_center_of_mass()
        self._solute_array -= Quantity(com, unit=nanometer)
        self.update_position()


    def export_solute_array(self, dir_name='', file_name='test.pdb'):
        '''
        アラニンジペプチドのpdbファイルを書き出す
        :param file_name:
        :return:
        '''
        coor = self._solute_array.in_units_of(angstrom)._value
        f = open(dir_name + file_name, 'w')
        f.write(f'MODEL        1\n\
HETATM    1  C   ACE A   1     {coor[0][0]:>6.3f}  {coor[0][1]:>6.3f}  {coor[0][2]:>6.3f}   1.00  0.00           C\n\
HETATM    2  O   ACE A   1     {coor[1][0]:>6.3f}  {coor[1][1]:>6.3f}  {coor[1][2]:>6.3f}   1.00  0.00           O\n\
HETATM    3  CH3 ACE A   1     {coor[2][0]:>6.3f}  {coor[2][1]:>6.3f}  {coor[2][2]:>6.3f}   1.00  0.00           C\n\
HETATM    4  H1  ACE A   1     {coor[3][0]:>6.3f}  {coor[3][1]:>6.3f}  {coor[3][2]:>6.3f}   1.00  0.00           H\n\
HETATM    5  H2  ACE A   1     {coor[4][0]:>6.3f}  {coor[4][1]:>6.3f}  {coor[4][2]:>6.3f}   1.00  0.00           H\n\
HETATM    6  H3  ACE A   1     {coor[5][0]:>6.3f}  {coor[5][1]:>6.3f}  {coor[5][2]:>6.3f}   1.00  0.00           H\n\
ATOM      7  N   ALA A   2     {coor[6][0]:>6.3f}  {coor[6][1]:>6.3f}  {coor[6][2]:>6.3f}   1.00  0.00           N\n\
ATOM      8  CA  ALA A   2     {coor[7][0]:>6.3f}  {coor[7][1]:>6.3f}  {coor[7][2]:>6.3f}   1.00  0.00           C\n\
ATOM      9  C   ALA A   2     {coor[8][0]:>6.3f}  {coor[8][1]:>6.3f}  {coor[8][2]:>6.3f}   1.00  0.00           C\n\
ATOM     10  O   ALA A   2     {coor[9][0]:>6.3f}  {coor[9][1]:>6.3f}  {coor[9][2]:>6.3f}   1.00  0.00           O\n\
ATOM     11  CB  ALA A   2     {coor[10][0]:>6.3f}  {coor[10][1]:>6.3f}  {coor[10][2]:>6.3f}   1.00  0.00           C\n\
ATOM     12  H   ALA A   2     {coor[11][0]:>6.3f}  {coor[11][1]:>6.3f}  {coor[11][2]:>6.3f}   1.00  0.00           H\n\
ATOM     13  HA  ALA A   2     {coor[12][0]:>6.3f}  {coor[12][1]:>6.3f}  {coor[12][2]:>6.3f}   1.00  0.00           H\n\
ATOM     14  HB1 ALA A   2     {coor[13][0]:>6.3f}  {coor[13][1]:>6.3f}  {coor[13][2]:>6.3f}   1.00  0.00           H\n\
ATOM     15  HB2 ALA A   2     {coor[14][0]:>6.3f}  {coor[14][1]:>6.3f}  {coor[14][2]:>6.3f}   1.00  0.00           H\n\
ATOM     16  HB3 ALA A   2     {coor[15][0]:>6.3f}  {coor[15][1]:>6.3f}  {coor[15][2]:>6.3f}   1.00  0.00           H\n\
HETATM   17  N   NME A   3     {coor[16][0]:>6.3f}  {coor[16][1]:>6.3f}  {coor[16][2]:>6.3f}   1.00  0.00           N\n\
HETATM   18  C   NME A   3     {coor[17][0]:>6.3f}  {coor[17][1]:>6.3f}  {coor[17][2]:>6.3f}   1.00  0.00           C\n\
HETATM   19  H   NME A   3     {coor[18][0]:>6.3f}  {coor[18][1]:>6.3f}  {coor[18][2]:>6.3f}   1.00  0.00           H\n\
HETATM   20  H1  NME A   3     {coor[19][0]:>6.3f}  {coor[19][1]:>6.3f}  {coor[19][2]:>6.3f}   1.00  0.00           H\n\
HETATM   21  H2  NME A   3     {coor[20][0]:>6.3f}  {coor[20][1]:>6.3f}  {coor[20][2]:>6.3f}   1.00  0.00           H\n\
HETATM   22  H3  NME A   3     {coor[21][0]:>6.3f}  {coor[21][1]:>6.3f}  {coor[21][2]:>6.3f}   1.00  0.00           H\n\
TER      23      NME A   3\n\
ENDMDL\n\
CONECT    1    3    2    7\n\
CONECT    2    1\n\
CONECT    3    1    4    5    6\n\
CONECT    4    3\n\
CONECT    5    3\n\
CONECT    6    3\n\
CONECT    7    1\n\
CONECT    9   17\n\
CONECT   17    9   18   19\n\
CONECT   18   20   21   22   17\n\
CONECT   19   17\n\
CONECT   20   18\n\
CONECT   21   18\n\
CONECT   22   18\n\
END')
        f.close()


if __name__ == '__main__':
    start_time = time.time()
    adp = ADP(dir_name='./C7eq/', file_name='C7eq.pdb')
    print(adp.solute_array)
    print('-----------------')
    adp.solute_array[0][1] += 1
    print(adp.solute_array)
    elapsed_time = time.time() - start_time
    print(f'elapsed time is {elapsed_time} [s]')