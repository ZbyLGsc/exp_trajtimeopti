# !/usr/bin/env python
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import time


def convertToFloat(dtl):
    dtl = dtl.split('=')[1]
    dtl = dtl.strip('[')
    dtl = dtl.strip(']\n')
    dtl = dtl.split(' ')
    dtl2 = []
    for i in range(len(dtl)):
        if not dtl[i] == '':
            dtl2.append(dtl[i])
    for i in range(len(dtl2)):
        dtl2[i] = float(dtl2[i])
    return dtl2


def main():
    # -----------------read data in file-----------------
    # get exp num
    num = input('Exp num?')
    ct_idx = num / 500
    ct_par = [[1, 2, 1], [1, 2, 2], [1, 2, 3], [1, 1, 1], [1, 3, 1], [2, 2, 1], [3, 2, 1]]
    file_name = 'data/exp_data'+str(ct_par[ct_idx][0]) + str(ct_par[ct_idx][1])+str(ct_par[ct_idx][2])+'.txt'
    # show trajectory screen shoot
    pic_name = 'pictures/'+str(num)+'.png'
    pic = mpimg.imread(pic_name)
    pic_crop = pic[190:-1, 40:-500, :]
    plt.axis('off')
    plt.imshow(pic_crop)

    # find exp data in data file
    df = file(file_name, 'r')
    # df = file('/home/zby/workspaces/exp_data.txt', 'r')
    while True:
        if df.readline().count('exp'+str(num)+'_start'):
            break
        else:
            continue
    # get data of traj 1
    v1t = convertToFloat(df.readline())
    v1x = convertToFloat(df.readline())
    v1y = convertToFloat(df.readline())
    v1z = convertToFloat(df.readline())
    a1t = convertToFloat(df.readline())
    a1x = convertToFloat(df.readline())
    a1y = convertToFloat(df.readline())
    a1z = convertToFloat(df.readline())
    c1t = convertToFloat(df.readline())
    c1 = convertToFloat(df.readline())
    print '\n-----Traj 1-----:'
    for i in range(10):
        print df.readline(),
    # get data of traj 2
    v2t = convertToFloat(df.readline())
    v2x = convertToFloat(df.readline())
    v2y = convertToFloat(df.readline())
    v2z = convertToFloat(df.readline())
    a2t = convertToFloat(df.readline())
    a2x = convertToFloat(df.readline())
    a2y = convertToFloat(df.readline())
    a2z = convertToFloat(df.readline())
    c2t = convertToFloat(df.readline())
    c2 = convertToFloat(df.readline())
    print '\n-----Traj 2-----:'
    for i in range(13):
        print df.readline(),
    df.close()

    # ---------------draw all data--------------------
    # create plot
    plt.figure(figsize=(30, 8), dpi=70)
    tvm = 20.0
    if v1t[len(v1t)-1] >= v2t[len(v2t)-1]:
        tvm = v1t[len(v1t)-1]
    else:
        tvm = v2t[len(v2t)-1]
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
    if a1t[len(a1t)-1] >= a2t[len(a2t)-1]:
        tam = a1t[len(a1t)-1]
    else:
        tam = a2t[len(a2t)-1]
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
    if c1t[len(c1t)-1] >= c2t[len(c2t)-1]:
        tcm = c1t[len(c1t)-1]
    else:
        tcm = c2t[len(c2t)-1]
    pc = plt.subplot(144)
    pc.axis([0.0, tcm, 0.0, 5.0])
    pc.set_xlabel("Iter time/ms", fontsize=14)
    pc.set_ylabel("Cost", fontsize=14)
    # draw x y z vel
    pvx.plot(v1t, v1x, 'r-', linewidth=2)
    pvx.plot(v2t, v2x, 'b-', linewidth=2)
    pvy.plot(v1t, v1y, 'r-', linewidth=2)
    pvy.plot(v2t, v2y, 'b-', linewidth=2)
    pvz.plot(v1t, v1z, 'r-', linewidth=2)
    pvz.plot(v2t, v2z, 'b-', linewidth=2)
    # draw vel bounding
    vb = []
    for i in range(len(v1t)):
        vb.append(2.0)
    pvx.plot(v1t, vb, 'g--', linewidth=2)
    pvy.plot(v1t, vb, 'g--', linewidth=2)
    pvz.plot(v1t, vb, 'g--', linewidth=2)
    vbmn = []
    for i in range(len(vb)):
        vbmn.append(-vb[i])
    pvx.plot(v1t, vbmn, 'g--', linewidth=2)
    pvy.plot(v1t, vbmn, 'g--', linewidth=2)
    pvz.plot(v1t, vbmn, 'g--', linewidth=2)
    # draw x y z acc
    pax.plot(a1t, a1x, 'r-', linewidth=2)
    pax.plot(a2t, a2x, 'b-', linewidth=2)
    pay.plot(a1t, a1y, 'r-', linewidth=2)
    pay.plot(a2t, a2y, 'b-', linewidth=2)
    paz.plot(a1t, a1z, 'r-', linewidth=2)
    paz.plot(a2t, a2z, 'b-', linewidth=2)
    # draw acc bounding
    ab = []
    for i in range(len(a1t)):
        ab.append(1.5)
    pax.plot(a1t, ab, 'g--', linewidth=2)
    pay.plot(a1t, ab, 'g--', linewidth=2)
    paz.plot(a1t, ab, 'g--', linewidth=2)
    abmn = []
    for i in range(len(ab)):
        abmn.append(-ab[i])
    pax.plot(a1t, abmn, 'g--', linewidth=2)
    pay.plot(a1t, abmn, 'g--', linewidth=2)
    paz.plot(a1t, abmn, 'g--', linewidth=2)
    # draw cost history
    pc.plot(c1t, c1, 'r-', linewidth=2)
    pc.plot(c2t, c2, 'b-', linewidth=2)
    plt.show()


if __name__ == '__main__':

    main()
