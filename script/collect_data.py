# !/usr/bin/env python
import rospy
from geometry_msgs.msg import PoseStamped
import matplotlib.pyplot as plt
import numpy as np
import pyscreenshot as pysc
import time

global vel1_data
global acc1_data
global cost_data
global len_data
global vel2_data
global acc2_data
vel1_data = []
acc1_data = []
cost_data = []
len_data = []
vel2_data = []
acc2_data = []
global can_draw
can_draw = False
global map_used_time
map_used_time = 0


def velCallback(msg):
    global vel1_data
    global can_draw
    # clear all data if receive all zero
    eps = 1e-3
    if (abs(msg.pose.position.x) < eps)and(abs(msg.pose.position.y) < eps)and(abs(msg.pose.position.z) < eps):
        vel1_data = []
        print 'receive new vel'
        can_draw = False
    # print and draw data if receive all one
    elif (abs(msg.pose.position.x - 1) < eps and
          abs(msg.pose.position.y - 1) < eps and abs(msg.pose.position.z - 1) < eps):
        print 'vel ok'
    else:
        # x as time, y and z for two domains' vel
        vel1_data.append(msg.pose.position.x)
        vel1_data.append(msg.pose.position.y)
        vel1_data.append(msg.pose.position.z)
        vel1_data.append(msg.pose.orientation.w)
        vel1_data.append(2.0)
        # print 'vel1:', len(vel1_data), vel1_data[len(vel1_data)-4:len(vel1_data)]


def accCallback(msg):
    global acc1_data
    # clear all data if receive all zero
    eps = 1e-3
    if (abs(msg.pose.position.x) < eps)and(abs(msg.pose.position.y) < eps)and(
            abs(msg.pose.position.z) < eps):
        acc1_data = []
        print 'receive new acc'
    # print and draw data if receive all one
    elif (abs(msg.pose.position.x - 1) < eps and
          abs(msg.pose.position.y - 1) < eps and abs(msg.pose.position.z - 1) < eps):
        print 'acc ok'
    else:
        # x as time, y and z for two domains' acc
        acc1_data.append(msg.pose.position.x)
        acc1_data.append(msg.pose.position.y)
        acc1_data.append(msg.pose.position.z)
        acc1_data.append(msg.pose.orientation.w)
        acc1_data.append(1.5)
        # print 'acc1:', len(acc1_data), acc1_data[len(acc1_data)-4:len(acc1_data)]


def vel2Callback(msg):
    global vel2_data
    # clear all data if receive all zero
    eps = 1e-3
    if (abs(msg.pose.position.x) < eps)and(abs(msg.pose.position.y) < eps)and(abs(msg.pose.position.z) < eps):
        vel2_data = []
        print 'receive new vel2'
    # print and draw data if receive all one
    elif (abs(msg.pose.position.x - 1) < eps and
          abs(msg.pose.position.y - 1) < eps and abs(msg.pose.position.z - 1) < eps):
        print 'vel2 ok'
    else:
        # x as time, y and z for two domains' vel
        vel2_data.append(msg.pose.position.x)
        vel2_data.append(msg.pose.position.y)
        vel2_data.append(msg.pose.position.z)
        vel2_data.append(msg.pose.orientation.w)
        # print 'vel2:', len(vel2_data), vel2_data[len(vel2_data)-4:len(vel2_data)]


def acc2Callback(msg):
    global acc2_data
    global can_draw
    # clear all data if receive all zero
    eps = 1e-3
    if (abs(msg.pose.position.x) < eps)and(abs(msg.pose.position.y) < eps)and(
            abs(msg.pose.position.z) < eps):
        acc2_data = []
        print 'receive new acc2'
    # print and draw data if receive all one
    elif (abs(msg.pose.position.x - 1) < eps and
          abs(msg.pose.position.y - 1) < eps and abs(msg.pose.position.z - 1) < eps):
        print 'acc2 ok'
        can_draw = True
    else:
        # x as time, y and z for two domains' acc
        acc2_data.append(msg.pose.position.x)
        acc2_data.append(msg.pose.position.y)
        acc2_data.append(msg.pose.position.z)
        acc2_data.append(msg.pose.orientation.w)
        # print 'acc2:', len(acc2_data), acc2_data[len(acc2_data)-4:len(acc2_data)]


def costCallback(msg):
    global cost_data
    # clear all data if receive all zero
    eps = 1e-5
    if (abs(msg.pose.position.x) < eps)and(abs(msg.pose.position.y) < eps)and(
            abs(msg.pose.position.z) < eps):
        cost_data = []
        print 'receive new cost'
    # print and draw data if receive all one
    elif ((abs(msg.pose.position.x - 1) < eps)and
          (abs(msg.pose.position.y - 1) < eps)and(abs(msg.pose.position.z - 1) < eps)):
        print 'cost ok'
    else:
        # x as time, y and z for two domains' cost
        cost_data.append(msg.pose.position.x)
        cost_data.append(msg.pose.position.y)
        cost_data.append(msg.pose.orientation.x)
        cost_data.append(msg.pose.orientation.y)
        # print 'cost:', len(cost_data), cost_data[len(cost_data)-4:len(cost_data)]


def lenCallback(msg):
    global len_data
    # clear all data if receive all zero
    eps = 1e-5
    if (abs(msg.pose.position.x) < eps)and(abs(msg.pose.position.y) < eps)and(
            abs(msg.pose.position.z) < eps):
        len_data = []
        print 'receive new len'
    # print and draw data if receive all one
    elif ((abs(msg.pose.position.x - 1) < eps)and
          (abs(msg.pose.position.y - 1) < eps)and(abs(msg.pose.position.z - 1) < eps)):
        print 'len ok'
    else:
        # x as time, y and z for two domains' cost
        len_data.append(msg.pose.position.x)
        len_data.append(msg.pose.position.y)
        len_data.append(msg.pose.position.z)
        # print 'len:', len_data


def main():
    print 'initial...'
    rospy.init_node("traj_test_drawer")
    print 'init ok...'
    # init subscriber
    rospy.Subscriber("time_optimal_traj/trajtest/vel", PoseStamped, velCallback)
    rospy.Subscriber("time_optimal_traj/trajtest/acc", PoseStamped, accCallback)
    rospy.Subscriber("time_optimal_traj/trajtest/vel2", PoseStamped, vel2Callback, queue_size=50)
    rospy.Subscriber("time_optimal_traj/trajtest/acc2", PoseStamped, acc2Callback, queue_size=50)
    rospy.Subscriber("time_optimal_traj/trajtest/cost", PoseStamped, costCallback)
    rospy.Subscriber("time_optimal_traj/trajtest/len", PoseStamped, lenCallback)
    print 'subscriber ready...'
    # init publisher
    nexp_pub = rospy.Publisher("time_optimal_traj/trajtest/nexp", PoseStamped, queue_size=1)
    nmap_pub = rospy.Publisher("/random_map/trajtest/nmap", PoseStamped, queue_size=1)
    time.sleep(0.5)
    print 'publisher ready...'
    # create subplot for vel, acc and cost
    plt.figure(figsize=(15, 6), dpi=50)
    plt.ion()
    plt.show()
    print 'matplotlib ready...'
    # ---------------main control vars-------------------------------
    global can_draw, vel1_data, acc1_data, cost_data, len_data
    can_draw = False
    root_dir = '/home/zby/Downloads/exp3/time_opti_exp1-2000/'
    # control important parameter
    weight_num = 1  # 1,2,3:0.5,1.5,4.5
    dense_num = 1  # 1,2,3:30,50,70
    grid_num = 1  # 1,2,3:30,20,10
    # weight + dense + grid, '141', '111', '151', '121', '131', '161'
    # ct_par = [[1, 2, 1], [1, 2, 2], [1, 2, 3], [1, 1, 1], [1, 3, 1], [2, 2, 1], [3, 2, 1],
    #           [1, 3, 1], [1, 2, 4], [1, 2, 5], [1, 2, 6], [4, 2, 1], [5, 2, 1], [6, 2, 1],
    #           [1, 4, 1], [1, 5, 1], [1, 6, 1], [1, 1, 1], [1, 2, 1], [1, 3, 1], [1, 1, 1],
    #           [1, 4, 1]]
    # ct_par = [[1, 4, 1], [1, 1, 1], [1, 5, 1], [1, 2, 1], [1, 3, 1], [1, 6, 1]]
    ct_par = [[1, 2, 1]]
    ct_idx = 0  # from 0 to 19
    file_name = root_dir+'data/exp_data'+str(ct_par[ct_idx][0]) + \
        str(ct_par[ct_idx][1])+str(ct_par[ct_idx][2])+'.txt'
    pic_dir = root_dir+'pictures/'
    exp_num = 0
    # start exp
    # ps = PoseStamped()
    # ps.pose.position.x = ct_par[ct_idx][1]
    # nmap_pub.publish(ps)
    # time.sleep(2.0)
    # ps.pose.position.x = ct_par[ct_idx][2]
    # ps.pose.position.y = ct_par[ct_idx][0]
    # nexp_pub.publish(ps)
    wait_time = 0
    print 'start exp!'
    print 'not able to draw yet...wait for new exp data...'
    while not rospy.is_shutdown():
        if can_draw:
            # -------------------sec1:collect experiment data and save to disk------------------------
            print 'can draw now, preparing...'
            # sort the array
            vt = []
            v1x = []
            v1y = []
            v1z = []
            vb = []
            at = []
            a1x = []
            a1y = []
            a1z = []
            ab = []
            ct1 = []
            ct2 = []
            c1 = []
            c2 = []
            vt2 = []
            v2x = []
            v2y = []
            v2z = []
            at2 = []
            a2x = []
            a2y = []
            a2z = []
            eps = 1e-4
            np.set_printoptions(linewidth='nan', threshold='nan')
            for i in range(len(vel1_data)):
                if abs(i % 5) <= eps:
                    vt.append(vel1_data[i])
                elif abs(i % 5-1) <= eps:
                    v1x.append(vel1_data[i])
                elif abs(i % 5-2) <= eps:
                    v1y.append(vel1_data[i])
                elif abs(i % 5-3) <= eps:
                    v1z.append(vel1_data[i])
                elif abs(i % 5-4) <= eps:
                    vb.append(vel1_data[i])
            vt = np.around(vt, 4)
            v1x = np.around(v1x, 4)
            v1y = np.around(v1y, 4)
            v1z = np.around(v1z, 4)
            for i in range(len(acc1_data)):
                if abs(i % 5) <= eps:
                    at.append(acc1_data[i])
                elif abs(i % 5-1) <= eps:
                    a1x.append(acc1_data[i])
                elif abs(i % 5-2) <= eps:
                    a1y.append(acc1_data[i])
                elif abs(i % 5-3) <= eps:
                    a1z.append(acc1_data[i])
                elif abs(i % 5-4) <= eps:
                    ab.append(acc1_data[i])
            at = np.around(at, 4)
            a1x = np.around(a1x, 4)
            a1y = np.around(a1y, 4)
            a1z = np.around(a1z, 4)
            for i in range(len(cost_data)):
                if abs(i % 4) <= eps:
                    ct1.append(cost_data[i])
                elif abs(i % 4-1) <= eps:
                    c1.append(cost_data[i])
                elif abs(i % 4-2) <= eps:
                    ct2.append(cost_data[i])
                elif abs(i % 4-3) <= eps:
                    c2.append(cost_data[i])
            ct1 = np.around(ct1, 8)
            c1 = np.around(c1, 8)
            ct2 = np.around(ct2, 8)
            c2 = np.around(c2, 8)
            for i in range(len(vel2_data)):
                if abs(i % 4) <= eps:
                    vt2.append(vel2_data[i])
                elif abs(i % 4-1) <= eps:
                    v2x.append(vel2_data[i])
                elif abs(i % 4-2) <= eps:
                    v2y.append(vel2_data[i])
                elif abs(i % 4-3) <= eps:
                    v2z.append(vel2_data[i])
            vt2 = np.around(vt2, 4)
            v2x = np.around(v2x, 4)
            v2y = np.around(v2y, 4)
            v2z = np.around(v2z, 4)
            for i in range(len(acc2_data)):
                if abs(i % 4) <= eps:
                    at2.append(acc2_data[i])
                elif abs(i % 4-1) <= eps:
                    a2x.append(acc2_data[i])
                elif abs(i % 4-2) <= eps:
                    a2y.append(acc2_data[i])
                elif abs(i % 4-3) <= eps:
                    a2z.append(acc2_data[i])
            at2 = np.around(at2, 4)
            a2x = np.around(a2x, 4)
            a2y = np.around(a2y, 4)
            a2z = np.around(a2z, 4)

            plt.clf()
            tvm = 20.0
            if vt[len(vt)-1] >= vt2[len(vt2)-1]:
                tvm = vt[len(vt)-1]
            else:
                tvm = vt2[len(vt2)-1]
            pvx = plt.subplot(241)
            pvx.axis([0.0, tvm, -3.5, 3.5])
            pvx.set_xlabel("Time/s", fontsize=14)
            pvx.set_ylabel("Vel x", fontsize=14)
            pvy = plt.subplot(242)
            pvy.axis([0.0, tvm, -3.5, 3.5])
            pvy.set_xlabel("Time/s", fontsize=14)
            pvy.set_ylabel("Vel y", fontsize=14)
            pvz = plt.subplot(243)
            pvz.axis([0.0, tvm, -3.5, 3.5])
            pvz.set_xlabel("Time/s", fontsize=14)
            pvz.set_ylabel("Vel z", fontsize=14)

            tam = 20.0
            if at[len(at)-1] >= at2[len(at2)-1]:
                tam = at[len(at)-1]
            else:
                tam = at2[len(at2)-1]
            pax = plt.subplot(245)
            pax.axis([0.0, tam, -3.0, 3.0])
            pax.set_xlabel("Time/s", fontsize=14)
            pax.set_ylabel("Acc x", fontsize=14)
            pay = plt.subplot(246)
            pay.axis([0.0, tam, -3.0, 3.0])
            pay.set_xlabel("Time/s", fontsize=14)
            pay.set_ylabel("Acc y", fontsize=14)
            paz = plt.subplot(247)
            paz.axis([0.0, tam, -3.0, 3.0])
            paz.set_xlabel("Time/s", fontsize=14)
            paz.set_ylabel("Acc z", fontsize=14)

            tcm = 200
            if ct1[len(ct1)-1] >= ct2[len(ct2)-1]:
                tcm = ct1[len(ct1)-1]
            else:
                tcm = ct2[len(ct2)-1]
            tcm = 1000 * tcm
            for i in range(len(ct1)):
                ct1[i] = 1000 * ct1[i]
            for i in range(len(ct2)):
                ct2[i] = 1000 * ct2[i]
            pc = plt.subplot(144)
            pc.axis([0.0, tcm, 0.0, 5.0])
            pc.set_xlabel("Iter time/ms", fontsize=14)
            pc.set_ylabel("Cost", fontsize=14)
            # draw x y z vel
            pvx.plot(vt, v1x, 'r-', linewidth=2)
            pvx.plot(vt2, v2x, 'b-', linewidth=2)
            pvy.plot(vt, v1y, 'r-', linewidth=2)
            pvy.plot(vt2, v2y, 'b-', linewidth=2)
            pvz.plot(vt, v1z, 'r-', linewidth=2)
            pvz.plot(vt2, v2z, 'b-', linewidth=2)
            # draw vel bounding
            pvx.plot(vt, vb, 'g--', linewidth=2)
            pvy.plot(vt, vb, 'g--', linewidth=2)
            pvz.plot(vt, vb, 'g--', linewidth=2)
            vbmn = []
            for i in range(len(vb)):
                vbmn.append(-vb[i])
            pvx.plot(vt, vbmn, 'g--', linewidth=2)
            pvy.plot(vt, vbmn, 'g--', linewidth=2)
            pvz.plot(vt, vbmn, 'g--', linewidth=2)
            # draw x y z acc
            pax.plot(at, a1x, 'r-', linewidth=2)
            pax.plot(at2, a2x, 'b-', linewidth=2)
            pay.plot(at, a1y, 'r-', linewidth=2)
            pay.plot(at2, a2y, 'b-', linewidth=2)
            paz.plot(at, a1z, 'r-', linewidth=2)
            paz.plot(at2, a2z, 'b-', linewidth=2)
            # draw acc bounding
            pax.plot(at, ab, 'g--', linewidth=2)
            pay.plot(at, ab, 'g--', linewidth=2)
            paz.plot(at, ab, 'g--', linewidth=2)
            abmn = []
            for i in range(len(ab)):
                abmn.append(-ab[i])
            pax.plot(at, abmn, 'g--', linewidth=2)
            pay.plot(at, abmn, 'g--', linewidth=2)
            paz.plot(at, abmn, 'g--', linewidth=2)
            # draw cost history
            pc.plot(ct1, c1, 'r-', linewidth=2)
            pc.plot(ct2, c2, 'b-', linewidth=2)
            plt.draw()
            plt.pause(2.0)

            # write data to file
            # wtf = raw_input('Save to file or not?[y/n]')
            wtf = 'y'
            if wtf == 'y':
                # num = raw_input('Exp num?')
                num = exp_num + 1
                # dtf = file('/home/zby/Downloads/exp3/time_opti_exp/exp_data.txt', 'a')
                dtf = file(file_name, 'a')
                dtf.write('exp'+str(num)+'_start----------------------------------\n')
                # write vel1, acc1, cost1
                print 'write traj 1...'
                dtf.write('v1t='+str(vt)+'\n')
                dtf.write('v1x='+str(v1x)+'\n')
                dtf.write('v1y='+str(v1y)+'\n')
                dtf.write('v1z='+str(v1z)+'\n')
                dtf.write('a1t='+str(at)+'\n')
                dtf.write('a1x='+str(a1x)+'\n')
                dtf.write('a1y='+str(a1y)+'\n')
                dtf.write('a1z='+str(a1z)+'\n')
                dtf.write('c1t='+str(ct1)+'\n')
                dtf.write('c1='+str(c1)+'\n')
                # write exceed ratio of traj 1
                v1x = np.max(np.abs(np.array(v1x)))
                ecr_v1x = np.max([v1x - 2.0, 0.0])/2.0
                v1y = np.max(np.abs(np.array(v1y)))
                ecr_v1y = np.max([v1y - 2.0, 0.0])/2.0
                v1z = np.max(np.abs(np.array(v1z)))
                ecr_v1z = np.max([v1z - 2.0, 0.0])/2.0
                a1x = np.max(np.abs(np.array(a1x)))
                ecr_a1x = np.max([a1x - 1.5, 0.0])/1.5
                a1y = np.max(np.abs(np.array(a1y)))
                ecr_a1y = np.max([a1y - 1.5, 0.0])/1.5
                a1z = np.max(np.abs(np.array(a1z)))
                ecr_a1z = np.max([a1z - 1.5, 0.0])/1.5
                dtf.write('excess ratio v1x:'+str(ecr_v1x)+'\n')
                dtf.write('excess ratio v1y:'+str(ecr_v1y)+'\n')
                dtf.write('excess ratio v1z:'+str(ecr_v1z)+'\n')
                dtf.write('excess ratio a1x:'+str(ecr_a1x)+'\n')
                dtf.write('excess ratio a1y:'+str(ecr_a1y)+'\n')
                dtf.write('excess ratio a1z:'+str(ecr_a1z)+'\n')
                # write length1, time1, meanv1
                dtf.write('len1:'+str(len_data[0])+'\n')
                dtf.write('time1:'+str(vt[len(vt)-1])+'\n')
                dtf.write('mean_v1:'+str(len_data[0]/vt[len(vt)-1])+'\n')
                # smoothness of traj 1:
                minc1 = np.min(c1)
                if minc1 > 1.0:
                    dtf.write('success1:'+str(0)+'\n')
                else:
                    dtf.write('success1:'+str(1)+'\n')

                # write vel2, acc2, cost2
                print 'write traj 2...'
                dtf.write('v2t='+str(vt2)+'\n')
                dtf.write('v2x='+str(v2x)+'\n')
                dtf.write('v2y='+str(v2y)+'\n')
                dtf.write('v2z='+str(v2z)+'\n')
                dtf.write('a2t='+str(at2)+'\n')
                dtf.write('a2x='+str(a2x)+'\n')
                dtf.write('a2y='+str(a2y)+'\n')
                dtf.write('a2z='+str(a2z)+'\n')
                dtf.write('c2t='+str(ct2)+'\n')
                dtf.write('c2='+str(c2)+'\n')
                # excess ratio of traj 2
                v2x = np.max(np.abs(np.array(v2x)))
                ecr_v2x = np.max([v2x - 2.0, 0.0])/2.0
                v2y = np.max(np.abs(np.array(v2y)))
                ecr_v2y = np.max([v2y - 2.0, 0.0])/2.0
                v2z = np.max(np.abs(np.array(v2z)))
                ecr_v2z = np.max([v2z - 2.0, 0.0])/2.0
                a2x = np.max(np.abs(np.array(a2x)))
                ecr_a2x = np.max([a2x - 1.5, 0.0])/1.5
                a2y = np.max(np.abs(np.array(a2y)))
                ecr_a2y = np.max([a2y - 1.5, 0.0])/1.5
                a2z = np.max(np.abs(np.array(a2z)))
                ecr_a2z = np.max([a2z - 1.5, 0.0])/1.5
                dtf.write('excess ratio v2x:'+str(ecr_v2x)+'\n')
                dtf.write('excess ratio v2y:'+str(ecr_v2y)+'\n')
                dtf.write('excess ratio v2z:'+str(ecr_v2z)+'\n')
                dtf.write('excess ratio a2x:'+str(ecr_a2x)+'\n')
                dtf.write('excess ratio a2y:'+str(ecr_a2y)+'\n')
                dtf.write('excess ratio a2z:'+str(ecr_a2z)+'\n')
                # write length2, time2, meanv2 and minimum time solver's time
                dtf.write('len2:'+str(len_data[1])+'\n')
                dtf.write('time2:'+str(vt2[len(vt2)-1])+'\n')
                dtf.write('mean_v2:'+str(len_data[1]/vt2[len(vt2)-1])+'\n')
                dtf.write('Min Solver Time:'+str(len_data[2])+'\n')
                dtf.write('Time per meter:'+str(len_data[2]/len_data[1])+'\n')
                # smoothness of traj 1:
                minc2 = np.min(c2)
                if minc2 > 1.0:
                    dtf.write('success2:'+str(0)+'\n')
                else:
                    dtf.write('success2:'+str(1)+'\n')
                # comparision of mean velocity
                dtf.write('Velocity ratio:' +
                          str((len_data[1]/vt2[len(vt2)-1])/(len_data[0]/vt[len(vt)-1]))+'\n')
                dtf.write('--------------------------------exp'+str(num)+'_end\n\n\n')
                dtf.close()
                print 'Write to file ok, take screen shot...'
                # screen shot of the trajectory
                img = pysc.grab()
                img.save(pic_dir+str(exp_num+1)+'.png')
                # y = raw_input('Continue?[Y/n]')
            can_draw = False
            # ---------------------sec2:prepare for new experiment----------------------
            # 500 exp/ group
            # exp_num += 1
            # exp_group_num = 4
            # print exp_num % exp_group_num
            # if abs(exp_num % exp_group_num) < 1e-3:
            #     print '500 exp are done!'
            #     ct_idx += 1
            #     if ct_idx > 1:
            #         print 'All exp done!!'
            #         return
            #     ps = PoseStamped()
            #     ps.pose.position.x = ct_par[ct_idx][1]
            #     nmap_pub.publish(ps)
            # file_name = root_dir+'data/exp_data'+str(ct_par[ct_idx][0]) + \
            #     str(ct_par[ct_idx][1])+str(ct_par[ct_idx][2])+'.txt'
            # print 'exp num:', exp_num
            # # update map after
            # # when time_optimal_traj restart, map need to be republish
            # ps = PoseStamped()
            # if abs(exp_num % 100) < 1e-3:
            #     print 'Main node restart, publish new map'
            #     time.sleep(2.0)
            #     ps.pose.position.x = ct_par[ct_idx][1]
            #     nmap_pub.publish(ps)
            #     time.sleep(3.0)
            # # when map is used for too many times, map need to be republish too.
            # global map_used_time
            # map_used_time += 1
            # if map_used_time > 9:
            #     time.sleep(2.0)
            #     # use different map density
            #     ps.pose.position.x = ct_par[ct_idx][1]
            #     nmap_pub.publish(ps)
            #     map_used_time = 0
            #     print 'update map'
            #     time.sleep(5.0)
            # # new experiment request, use different weight and different grid_num
            # ps.pose.position.x = ct_par[ct_idx][2]
            # ps.pose.position.y = ct_par[ct_idx][0]
            # nexp_pub.publish(ps)
            # wait_time = 0
            print '---------------start new exp!------------------'
            print 'not able to draw yet...wait for new exp data...'
        else:
            time.sleep(1.0)
            # if wait for too long, restart a new exp
            # wait_time += 1
            # print 'waiting time:', wait_time
            # if abs(wait_time % 5) <1e-3:
            #     ps = PoseStamped()
            #     ps.pose.position.x = ct_par[ct_idx][2]
            #     ps.pose.position.y = ct_par[ct_idx][0]
            #     nexp_pub.publish(ps)
            # if wait_time >= 150:
            #     print 'Main node may die unexpectedly, publish new map'
            #     time.sleep(5.0)
            #     ps.pose.position.x = ct_par[ct_idx][1]
            #     nmap_pub.publish(ps)
            #     time.sleep(5.0)
            #     wait_time = 0
                

    rospy.spin()


if __name__ == '__main__':

    try:
        main()
    except rospy.ROSInterruptException:
        pass
