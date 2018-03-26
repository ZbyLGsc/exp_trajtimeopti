# !/usr/bin/env python
import matplotlib.pyplot as plt
from draw_data import convertToFloat


def addFileEnd():
    f_idxs = ['111', '121', '131', '141', '151', '161']
    f_prefix = 'data/exp_data'
    f_suffix = '.txt'
    for id in f_idxs:
        f = file(f_prefix+id+f_suffix, 'a')
        f.write('\n')
        f.write('file_end')
        f.close()


def find50MSIndex(line):
    # convert to float list
    line = convertToFloat(line)
    # find the first time larger than 50.0
    idx = -1
    for i in range(len(line)):
        if line[i] > 50.0:
            idx = i
            break
    return idx


def getOptimizationCost(line, idx):
    line = convertToFloat(line)
    return line[idx] / 4.0 * 100.0


def compareGridNum():
    print 'compare different grid num...'
    f_idxs = ['126', '123', '125', '122', '124']
    f_prefix = 'data/exp_data'
    f_suffix = '.txt'
    # average of solving time and velocity ratio
    slv_time_list = []
    traj_len_list = []
    vel_ratio_list = []
    # traversing each file
    for idx in f_idxs:
        f_name = f_prefix + idx + f_suffix
        dtf = file(f_name, 'r')
        print 'open file: ', f_name
        slv_time = []
        traj_len = []
        vel_ratio = []
        while True:
            line = dtf.readline()
            # stop traversing if reach end of file
            if line.count('file_end'):
                print 'end of file: ', f_name
                break
            # average solving time and velocity ratio
            elif line.count('Min Solver Time'):
                slv_time.append(1000.0 * float(line.split(':')[1].strip('\n')))
            elif line.count('len2'):
                traj_len.append(float(line.split(':')[1].strip('\n')))
            elif line.count('Velocity ratio'):
                vel_ratio.append(float(line.split(':')[1].strip('\n')))
        # calculate average
        slv_time = sum(slv_time) / len(slv_time)
        traj_len = sum(traj_len) / len(traj_len)
        vel_ratio = sum(vel_ratio) / len(vel_ratio)

        slv_time_list.append(slv_time)
        traj_len_list.append(traj_len)
        vel_ratio_list.append(vel_ratio)
    # plot data
    print 'avg traj len:', traj_len_list
    gn = [5, 10, 15, 20, 25]
    plt.figure('solve time')
    plt.grid()
    plt.axis([5, 30, 0, 100])
    plt.xlabel('Grid num')
    plt.ylabel('Solve time/ms')
    plt.plot(gn, slv_time_list)

    plt.figure('vel ratio')
    plt.grid()
    plt.axis([5, 30, 0, 2])
    plt.xlabel('Grid num')
    plt.ylabel('vel ratio')
    plt.plot(gn, vel_ratio_list)
    plt.show()


def compareWeight():
    print 'compare different weight...'
    # 142563
    f_idxs = ['421', '221', '521', '621', '321']
    f_prefix = 'data/exp_data'
    f_suffix = '.txt'
    # average of velocity ratio
    traj_len_list = []
    vel_ratio_list = []
    # traversing each file
    for idx in f_idxs:
        f_name = f_prefix + idx + f_suffix
        dtf = file(f_name, 'r')
        print 'open file: ', f_name
        traj_len = []
        vel_ratio = []
        while True:
            line = dtf.readline()
            # stop traversing if reach end of file
            if line.count('file_end'):
                print 'end of file: ', f_name
                break
            # average velocity ratio
            elif line.count('len2'):
                traj_len.append(float(line.split(':')[1].strip('\n')))
            elif line.count('Velocity ratio'):
                vel_ratio.append(float(line.split(':')[1].strip('\n')))
        # calculate average
        traj_len = sum(traj_len) / len(traj_len)
        vel_ratio = sum(vel_ratio) / len(vel_ratio)

        traj_len_list.append(traj_len)
        vel_ratio_list.append(vel_ratio)
    # plot data
    print 'avg traj len:', traj_len_list
    we = [3.0, 5.0, 8.0, 11.0, 15.0]
    plt.figure('vel ratio')
    plt.grid()
    plt.axis([0, 15, 0, 2])
    plt.xlabel('Weight')
    plt.ylabel('vel ratio')
    plt.plot(we, vel_ratio_list)
    plt.show()


def compareDensity():
    print 'compare defferent density...'
    # specify data file
    f_idxs = ['141', '111', '151', '121', '131', '161']
    # f_idxs = ['111', '151', '121', '131', '161']
    f_prefix = 'data/exp_data'
    f_suffix = '.txt'
    # data to be collect
    cost1s_list = []
    cost2s_list = []
    exc_v1x_list = []
    exc_v1y_list = []
    exc_v1z_list = []
    exc_a1x_list = []
    exc_a1y_list = []
    exc_a1z_list = []
    exc_v2x_list = []
    exc_v2y_list = []
    exc_v2z_list = []
    exc_a2x_list = []
    exc_a2y_list = []
    exc_a2z_list = []
    vel_ratio_list = []
    vel1_list = []
    vel2_list = [] 
    exc_list = []
    # calculate 50ms cost and excess ratio of each 500 exp
    for idx in f_idxs:
        f_name = f_prefix + idx + f_suffix
        dtf = file(f_name, 'r')
        print 'open file: ', f_name
        # cost1, cost2, excess1, excess2
        cost1s = []
        cost2s = []
        exc_v1x = []
        exc_v1y = []
        exc_v1z = []
        exc_a1x = []
        exc_a1y = []
        exc_a1z = []
        exc_v2x = []
        exc_v2y = []
        exc_v2z = []
        exc_a2x = []
        exc_a2y = []
        exc_a2z = []
        vel_ratio = []
        vel1 = []
        vel2 = []
        # traverse the file
        while True:
            line = dtf.readline()
            # stop traversing if reach end of file
            if line.count('file_end'):
                print 'end of file: ', f_name
                break
            elif line.count('c1t='):
                idx_50 = find50MSIndex(line)
                c1_line = dtf.readline()
                cost1 = getOptimizationCost(c1_line, idx_50)
                cost1s.append(cost1)
            elif line.count('c2t='):
                idx_50 = find50MSIndex(line)
                c2_line = dtf.readline()
                cost2 = getOptimizationCost(c2_line, idx_50)
                cost2s.append(cost2)
            elif line.count('Velocity ratio'):
                vel_ratio.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio v1x'):
                exc_v1x.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio v1y'):
                exc_v1y.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio v1z'):
                exc_v1z.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio a1x'):
                exc_a1x.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio a1y'):
                exc_a1y.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio a1z'):
                exc_a1z.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio v2x'):
                exc_v2x.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio v2y'):
                exc_v2y.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio v2z'):
                exc_v2z.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio a2x'):
                exc_a2x.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio a2y'):
                exc_a2y.append(float(line.split(':')[1].strip('\n')))
            elif line.count('excess ratio a2z'):
                exc_a2z.append(float(line.split(':')[1].strip('\n')))
            elif line.count('mean_v1'):
                vel1.append(float(line.split(':')[1].strip('\n')))
            elif line.count('mean_v2'):
                vel2.append(float(line.split(':')[1].strip('\n')))
        # record excess num
        exn_v = 0
        exn_a = 0
        for i in range(len(exc_v1x)):
            if exc_v1x[i] > 0.0 or exc_v1y[i] > 0.0 or exc_a1x[i] > 0.0 or exc_a1y[i] > 0.0:
                exn_v += 1
        exc_list.append(exn_v / 5.0)
        # calculate average
        print 'calc average...'
        cost1_avrg = 100 - sum(cost1s)/len(cost1s)
        cost2_avrg = 100 - sum(cost2s)/len(cost2s)
        ev1x_avrg = sum(exc_v1x)/len(exc_v1x)
        ev1y_avrg = sum(exc_v1y)/len(exc_v1y)
        ev1z_avrg = sum(exc_v1z)/len(exc_v1z)
        ea1x_avrg = sum(exc_a1x)/len(exc_a1x)
        ea1y_avrg = sum(exc_a1y)/len(exc_a1y)
        ea1z_avrg = sum(exc_a1z)/len(exc_a1z)
        ev2x_avrg = sum(exc_v2x)/len(exc_v2x)
        ev2y_avrg = sum(exc_v2y)/len(exc_v2y)
        ev2z_avrg = sum(exc_v2z)/len(exc_v2z)
        ea2x_avrg = sum(exc_a2x)/len(exc_a2x)
        ea2y_avrg = sum(exc_a2y)/len(exc_a2y)
        ea2z_avrg = sum(exc_a2z)/len(exc_a2z)
        vel_avrg = sum(vel_ratio)/len(vel_ratio)
        vel1_avrg = sum(vel1)/len(vel1)
        vel2_avrg = sum(vel2)/len(vel2)
        # add to list
        cost1s_list.append(cost1_avrg)
        cost2s_list.append(cost2_avrg)
        exc_v1x_list.append(ev1x_avrg)
        exc_v1y_list.append(ev1y_avrg)
        exc_v1z_list.append(ev1z_avrg)
        exc_a1x_list.append(ea1x_avrg)
        exc_a1y_list.append(ea1y_avrg)
        exc_a1z_list.append(ea1z_avrg)
        exc_v2x_list.append(ev2x_avrg)
        exc_v2y_list.append(ev2y_avrg)
        exc_v2z_list.append(ev2z_avrg)
        exc_a2x_list.append(ea2x_avrg)
        exc_a2y_list.append(ea2y_avrg)
        exc_a2z_list.append(ea2z_avrg)
        vel_ratio_list.append(vel_avrg)
        vel1_list.append(vel1_avrg)
        vel2_list.append(vel2_avrg)
        dtf.close()
        print 'file: ', f_name, ' finish'
    # plot cost
    den = [20, 30, 40, 50, 60, 70]
    plt.figure('cost')
    plt.grid()
    plt.axis([20, 70, 0, 100])
    plt.ylabel("Cost reduction(%)", fontsize=14)
    plt.xlabel("Obstacle density(Trees/map)", fontsize=14)
    plt.plot(den, cost1s_list, 'r', linewidth=2, marker='s', label='Our previous method')
    plt.plot(den, cost2s_list, 'b', linewidth=2, marker='s', label='Our proposed method')
    # plot mean vel
    plt.twinx()
    plt.ylabel('Average Velocity(m/s)', fontsize=14)
    plt.ylim(0.5, 2.0)
    plt.plot(den, vel1_list, 'r', linewidth=2, linestyle='--')
    plt.plot(den, vel2_list, 'b', linewidth=2, linestyle='--')
    print 'exc:', exc_list
    # plt.savefig()
    # print 'den=', den
    # print 'cost1=', cost1s_list
    # print 'cost2=', cost2s_list
    # print 'vel1=', vel1_list
    # print 'vel2=', vel2_list
    # print 'vel_ratio=', vel_ratio_list
    # print 'vx=', exc_v1x_list
    # print 'vy=', exc_v1y_list
    # print 'ax=', exc_a1x_list
    # print 'ay=', exc_a1y_list
    # # plot excess
    # plt.figure('x vel')
    # plt.grid()
    # plt.axis([20, 70, 0, 0.5])
    # plt.xlabel('density')
    # plt.ylabel('excess vel x')
    # print 'v1x:', exc_v1x_list
    # plt.plot(den, exc_v1x_list, 'r', linewidth=2)
    # plt.plot(den, exc_v2x_list, 'b', linewidth=2)

    # plt.figure('y vel')
    # plt.grid()
    # plt.axis([20, 70, 0, 0.5])
    # plt.xlabel('density')
    # plt.ylabel('excess vel y')
    # plt.plot(den, exc_v1y_list, 'r', linewidth=2)
    # plt.plot(den, exc_v2y_list, 'b', linewidth=2)

    # plt.figure('x acc')
    # plt.grid()
    # plt.axis([20, 70, 0, 1.0])
    # plt.xlabel('density')
    # plt.ylabel('excess acc x')
    # plt.plot(den, exc_a1x_list, 'r', linewidth=2)
    # plt.plot(den, exc_a2x_list, 'b', linewidth=2)

    # plt.figure('y acc')
    # plt.grid()
    # plt.axis([20, 70, 0, 1.0])
    # plt.xlabel('density')
    # plt.ylabel('excess acc y')
    # plt.plot(den, exc_a1y_list, 'r', linewidth=2)
    # plt.plot(den, exc_a2y_list, 'b', linewidth=2)
    plt.show()


if __name__ == '__main__':
    # addFileEnd()
    compareDensity()
    # compareGridNum()
    # compareWeight()
