from merapy.graphics import Graphics

G_2_2 = {}
if 1:
    # I: [('O', 2), ('O', 4)]
    # O: [('Vp_1', 1), ('I', 1), ('V_1', 6), ('I', 2)]
    # V_1: [('Vp_1', 2), ('oo_2_3', 3), ('oo_2_3', 4), ('Vp_1', 5), ('Vp_1', 6), ('O', 3)]
    # Vp_1: [('O', 1), ('V_1', 1), ('oo_2_3', 1), ('oo_2_3', 2), ('V_1', 4), ('V_1', 5)]
    # oo_2_3: [('Vp_1', 3), ('Vp_1', 4), ('V_1', 2), ('V_1', 3)]
    G_2_2[0]=Graphics()
    G_2_2[0].graph_name='G_2_2[0]'
    G_2_2[0].size=5
    G_2_2[0].nodes[:]=-1
    G_2_2[0].edges[:, :]=-1
    G_2_2[0].names[0]="I"  #I
    G_2_2[0].nodes[0]=2
    G_2_2[0].edges[0:2,0] = [0,1]

    G_2_2[0].names[1]="oo"  #oo_2_3
    G_2_2[0].nodes[1]=4
    G_2_2[0].edges[0:4,1] = [2,3,4,5]

    G_2_2[0].names[2]="OO"  #O
    G_2_2[0].nodes[2]=4
    G_2_2[0].edges[0:4,2] = [6,0,7,1]

    G_2_2[0].names[3]="V"  #V_1
    G_2_2[0].nodes[3]=6
    G_2_2[0].edges[0:6,3] = [8,4,5,9,10,7]

    G_2_2[0].names[4]="Vp"  #Vp_1
    G_2_2[0].nodes[4]=6
    G_2_2[0].edges[0:6,4] = [6,8,2,3,9,10]

    # 6/8       2.02000002E+16
    G_2_2[0].contract_order=[0,2,3,4,1]  

    # I: [('O', 2), ('O', 4)]
    # O: [('Vp_1', 1), ('I', 1), ('V_1', 6), ('I', 2)]
    # V_1: [('Vp_1', 2), ('Vp_1', 3), ('oo_3_4', 3), ('oo_3_4', 4), ('Vp_1', 6), ('O', 3)]
    # Vp_1: [('O', 1), ('V_1', 1), ('V_1', 2), ('oo_3_4', 1), ('oo_3_4', 2), ('V_1', 5)]
    # oo_3_4: [('Vp_1', 4), ('Vp_1', 5), ('V_1', 3), ('V_1', 4)]
    G_2_2[1]=Graphics()
    G_2_2[1].graph_name='G_2_2[1]'
    G_2_2[1].size=5
    G_2_2[1].nodes[:]=-1
    G_2_2[1].edges[:, :]=-1
    G_2_2[1].names[0]="I"  #I
    G_2_2[1].nodes[0]=2
    G_2_2[1].edges[0:2,0] = [0,1]

    G_2_2[1].names[1]="OO"  #O
    G_2_2[1].nodes[1]=4
    G_2_2[1].edges[0:4,1] = [2,0,3,1]

    G_2_2[1].names[2]="V"  #V_1
    G_2_2[1].nodes[2]=6
    G_2_2[1].edges[0:6,2] = [4,5,6,7,8,3]

    G_2_2[1].names[3]="Vp"  #Vp_1
    G_2_2[1].nodes[3]=6
    G_2_2[1].edges[0:6,3] = [2,4,5,9,10,8]

    G_2_2[1].names[4]="oo"  #oo_3_4
    G_2_2[1].nodes[4]=4
    G_2_2[1].edges[0:4,4] = [9,10,6,7]

    # 6/8       4.04000002E+16
    G_2_2[1].contract_order=[0,1,2,3,4]  

    # OO: [('Vp_1', 1), ('Vp_2', 1), ('V_1', 6), ('V_2', 6)]
    # U_2: [('oo_4_5', 4), ('Up_2', 4), ('V_1', 5), ('V_2', 1)]
    # Up_2: [('Vp_1', 6), ('Vp_2', 2), ('oo_4_5', 2), ('U_2', 2)]
    # V_1: [('Vp_1', 2), ('Vp_1', 3), ('Vp_1', 4), ('oo_4_5', 3), ('U_2', 3), ('OO', 3)]
    # V_2: [('U_2', 4), ('Vp_2', 3), ('Vp_2', 4), ('Vp_2', 5), ('Vp_2', 6), ('OO', 4)]
    # Vp_1: [('OO', 1), ('V_1', 1), ('V_1', 2), ('V_1', 3), ('oo_4_5', 1), ('Up_2', 1)]
    # Vp_2: [('OO', 2), ('Up_2', 2), ('V_2', 2), ('V_2', 3), ('V_2', 4), ('V_2', 5)]
    # oo_4_5: [('Vp_1', 5), ('Up_2', 3), ('V_1', 4), ('U_2', 1)]
    G_2_2[2]=Graphics()
    G_2_2[2].graph_name='G_2_2[2]'
    G_2_2[2].size=8
    G_2_2[2].nodes[:]=-1
    G_2_2[2].edges[:, :]=-1
    G_2_2[2].names[0]="Up"  #Up_2
    G_2_2[2].nodes[0]=4
    G_2_2[2].edges[0:4,0] = [0,1,2,3]

    G_2_2[2].names[1]="OO"  #OO
    G_2_2[2].nodes[1]=4
    G_2_2[2].edges[0:4,1] = [4,5,6,7]

    G_2_2[2].names[2]="V"  #V_1
    G_2_2[2].nodes[2]=6
    G_2_2[2].edges[0:6,2] = [8,9,10,11,12,6]

    G_2_2[2].names[3]="V"  #V_2
    G_2_2[2].nodes[3]=6
    G_2_2[2].edges[0:6,3] = [13,14,15,16,17,7]

    G_2_2[2].names[4]="Vp"  #Vp_2
    G_2_2[2].nodes[4]=6
    G_2_2[2].edges[0:6,4] = [5,1,14,15,16,17]

    G_2_2[2].names[5]="Vp"  #Vp_1
    G_2_2[2].nodes[5]=6
    G_2_2[2].edges[0:6,5] = [4,8,9,10,18,0]

    G_2_2[2].names[6]="U"  #U_2
    G_2_2[2].nodes[6]=4
    G_2_2[2].edges[0:4,6] = [19,3,12,13]

    G_2_2[2].names[7]="oo"  #oo_4_5
    G_2_2[2].nodes[7]=4
    G_2_2[2].edges[0:4,7] = [18,2,11,19]

    # 6/9         4.070202E+18
    G_2_2[2].contract_order=[2,5,1,7,0,6,3,4]  

    # OO: [('Vp_1', 1), ('Vp_2', 1), ('V_1', 6), ('V_2', 6)]
    # U_2: [('oo_5_6', 3), ('oo_5_6', 4), ('V_1', 5), ('V_2', 1)]
    # Up_2: [('Vp_1', 6), ('Vp_2', 2), ('oo_5_6', 1), ('oo_5_6', 2)]
    # V_1: [('Vp_1', 2), ('Vp_1', 3), ('Vp_1', 4), ('Vp_1', 5), ('U_2', 3), ('OO', 3)]
    # V_2: [('U_2', 4), ('Vp_2', 3), ('Vp_2', 4), ('Vp_2', 5), ('Vp_2', 6), ('OO', 4)]
    # Vp_1: [('OO', 1), ('V_1', 1), ('V_1', 2), ('V_1', 3), ('V_1', 4), ('Up_2', 1)]
    # Vp_2: [('OO', 2), ('Up_2', 2), ('V_2', 2), ('V_2', 3), ('V_2', 4), ('V_2', 5)]
    # oo_5_6: [('Up_2', 3), ('Up_2', 4), ('U_2', 1), ('U_2', 2)]
    G_2_2[3]=Graphics()
    G_2_2[3].graph_name='G_2_2[3]'
    G_2_2[3].size=8
    G_2_2[3].nodes[:]=-1
    G_2_2[3].edges[:, :]=-1
    G_2_2[3].names[0]="Up"  #Up_2
    G_2_2[3].nodes[0]=4
    G_2_2[3].edges[0:4,0] = [0,1,2,3]

    G_2_2[3].names[1]="OO"  #OO
    G_2_2[3].nodes[1]=4
    G_2_2[3].edges[0:4,1] = [4,5,6,7]

    G_2_2[3].names[2]="V"  #V_1
    G_2_2[3].nodes[2]=6
    G_2_2[3].edges[0:6,2] = [8,9,10,11,12,6]

    G_2_2[3].names[3]="V"  #V_2
    G_2_2[3].nodes[3]=6
    G_2_2[3].edges[0:6,3] = [13,14,15,16,17,7]

    G_2_2[3].names[4]="Vp"  #Vp_2
    G_2_2[3].nodes[4]=6
    G_2_2[3].edges[0:6,4] = [5,1,14,15,16,17]

    G_2_2[3].names[5]="Vp"  #Vp_1
    G_2_2[3].nodes[5]=6
    G_2_2[3].edges[0:6,5] = [4,8,9,10,11,0]

    G_2_2[3].names[6]="U"  #U_2
    G_2_2[3].nodes[6]=4
    G_2_2[3].edges[0:4,6] = [18,19,12,13]

    G_2_2[3].names[7]="oo"  #oo_5_6
    G_2_2[3].nodes[7]=4
    G_2_2[3].edges[0:4,7] = [2,3,18,19]

    # 6/8           9.0303E+16
    G_2_2[3].contract_order=[2,5,1,0,7,6,3,4]  

    # OO: [('Vp_1', 1), ('Vp_2', 1), ('V_1', 6), ('V_2', 6)]
    # U_2: [('Up_2', 3), ('oo_6_7', 3), ('V_1', 5), ('V_2', 1)]
    # Up_2: [('Vp_1', 6), ('Vp_2', 2), ('U_2', 1), ('oo_6_7', 1)]
    # V_1: [('Vp_1', 2), ('Vp_1', 3), ('Vp_1', 4), ('Vp_1', 5), ('U_2', 3), ('OO', 3)]
    # V_2: [('U_2', 4), ('oo_6_7', 4), ('Vp_2', 4), ('Vp_2', 5), ('Vp_2', 6), ('OO', 4)]
    # Vp_1: [('OO', 1), ('V_1', 1), ('V_1', 2), ('V_1', 3), ('V_1', 4), ('Up_2', 1)]
    # Vp_2: [('OO', 2), ('Up_2', 2), ('oo_6_7', 2), ('V_2', 3), ('V_2', 4), ('V_2', 5)]
    # oo_6_7: [('Up_2', 4), ('Vp_2', 3), ('U_2', 2), ('V_2', 2)]
    G_2_2[4]=Graphics()
    G_2_2[4].graph_name='G_2_2[4]'
    G_2_2[4].size=8
    G_2_2[4].nodes[:]=-1
    G_2_2[4].edges[:, :]=-1
    G_2_2[4].names[0]="Up"  #Up_2
    G_2_2[4].nodes[0]=4
    G_2_2[4].edges[0:4,0] = [0,1,2,3]

    G_2_2[4].names[1]="OO"  #OO
    G_2_2[4].nodes[1]=4
    G_2_2[4].edges[0:4,1] = [4,5,6,7]

    G_2_2[4].names[2]="V"  #V_1
    G_2_2[4].nodes[2]=6
    G_2_2[4].edges[0:6,2] = [8,9,10,11,12,6]

    G_2_2[4].names[3]="V"  #V_2
    G_2_2[4].nodes[3]=6
    G_2_2[4].edges[0:6,3] = [13,14,15,16,17,7]

    G_2_2[4].names[4]="oo"  #oo_6_7
    G_2_2[4].nodes[4]=4
    G_2_2[4].edges[0:4,4] = [3,18,19,14]

    G_2_2[4].names[5]="Vp"  #Vp_2
    G_2_2[4].nodes[5]=6
    G_2_2[4].edges[0:6,5] = [5,1,18,15,16,17]

    G_2_2[4].names[6]="Vp"  #Vp_1
    G_2_2[4].nodes[6]=6
    G_2_2[4].edges[0:6,6] = [4,8,9,10,11,0]

    G_2_2[4].names[7]="U"  #U_2
    G_2_2[4].nodes[7]=4
    G_2_2[4].edges[0:4,7] = [2,19,12,13]

    # 6/9         4.100102E+18
    G_2_2[4].contract_order=[2,6,1,0,7,4,3,5]  

for i in G_2_2:
    G_2_2[i].weight = 1./5

