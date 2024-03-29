from merapy.graphics import Graphics

"""
gtop is the same as ternary, so wont define it again
"""

__all__ = ["G_2_2", "G_3_3", "G_3_2", "G_22_2", "G_22_3"]
G = {}

G_2_2 = {}
if 1:
    # OO: [('Vp_0', 1), ('Vp_1', 1), ('V_0', 3), ('V_1', 3)]
    # U_0: [('oo_0_1', 4), ('Up_0', 4), ('V_0', 2), ('V_1', 1)]
    # Up_0: [('Vp_0', 3), ('Vp_1', 2), ('oo_0_1', 2), ('U_0', 2)]
    # V_0: [('oo_0_1', 3), ('U_0', 3), ('OO', 3)]
    # V_1: [('U_0', 4), ('Vp_1', 3), ('OO', 4)]
    # Vp_0: [('OO', 1), ('oo_0_1', 1), ('Up_0', 1)]
    # Vp_1: [('OO', 2), ('Up_0', 2), ('V_1', 2)]
    # oo_0_1: [('Vp_0', 2), ('Up_0', 3), ('V_0', 1), ('U_0', 1)]
    G_2_2[0]=Graphics()
    G_2_2[0].graph_name='G_2_2[0]'
    G_2_2[0].size=8
    G_2_2[0].nodes[:]=-1
    G_2_2[0].edges[:, :]=-1
    G_2_2[0].names[0]="OO"  #OO
    G_2_2[0].nodes[0]=4
    G_2_2[0].edges[0:4,0] = [0,1,2,3]

    G_2_2[0].names[1]="Up"  #Up_0
    G_2_2[0].nodes[1]=4
    G_2_2[0].edges[0:4,1] = [4,5,6,7]

    G_2_2[0].names[2]="V"  #V_0
    G_2_2[0].nodes[2]=3
    G_2_2[0].edges[0:3,2] = [8,9,2]

    G_2_2[0].names[3]="V"  #V_1
    G_2_2[0].nodes[3]=3
    G_2_2[0].edges[0:3,3] = [10,11,3]

    G_2_2[0].names[4]="Vp"  #Vp_1
    G_2_2[0].nodes[4]=3
    G_2_2[0].edges[0:3,4] = [1,5,11]

    G_2_2[0].names[5]="Vp"  #Vp_0
    G_2_2[0].nodes[5]=3
    G_2_2[0].edges[0:3,5] = [0,12,4]

    G_2_2[0].names[6]="U"  #U_0
    G_2_2[0].nodes[6]=4
    G_2_2[0].edges[0:4,6] = [13,7,9,10]

    G_2_2[0].names[7]="oo"  #oo_0_1
    G_2_2[0].nodes[7]=4
    G_2_2[0].edges[0:4,7] = [12,6,8,13]

    # 6/7      1.010000002E+15
    G_2_2[0].weight = 0.250000000000000  
    G_2_2[0].contract_order=[5,7,1,2,6,3,0,4]  

    # OO: [('Vp_0', 1), ('Vp_1', 1), ('V_0', 3), ('V_1', 3)]
    # U_0: [('oo_1_2', 3), ('oo_1_2', 4), ('V_0', 2), ('V_1', 1)]
    # Up_0: [('Vp_0', 3), ('Vp_1', 2), ('oo_1_2', 1), ('oo_1_2', 2)]
    # V_0: [('Vp_0', 2), ('U_0', 3), ('OO', 3)]
    # V_1: [('U_0', 4), ('Vp_1', 3), ('OO', 4)]
    # Vp_0: [('OO', 1), ('V_0', 1), ('Up_0', 1)]
    # Vp_1: [('OO', 2), ('Up_0', 2), ('V_1', 2)]
    # oo_1_2: [('Up_0', 3), ('Up_0', 4), ('U_0', 1), ('U_0', 2)]
    G_2_2[1]=Graphics()
    G_2_2[1].graph_name='G_2_2[1]'
    G_2_2[1].size=8
    G_2_2[1].nodes[:]=-1
    G_2_2[1].edges[:, :]=-1
    G_2_2[1].names[0]="OO"  #OO
    G_2_2[1].nodes[0]=4
    G_2_2[1].edges[0:4,0] = [0,1,2,3]

    G_2_2[1].names[1]="Up"  #Up_0
    G_2_2[1].nodes[1]=4
    G_2_2[1].edges[0:4,1] = [4,5,6,7]

    G_2_2[1].names[2]="V"  #V_0
    G_2_2[1].nodes[2]=3
    G_2_2[1].edges[0:3,2] = [8,9,2]

    G_2_2[1].names[3]="V"  #V_1
    G_2_2[1].nodes[3]=3
    G_2_2[1].edges[0:3,3] = [10,11,3]

    G_2_2[1].names[4]="oo"  #oo_1_2
    G_2_2[1].nodes[4]=4
    G_2_2[1].edges[0:4,4] = [6,7,12,13]

    G_2_2[1].names[5]="Vp"  #Vp_1
    G_2_2[1].nodes[5]=3
    G_2_2[1].edges[0:3,5] = [1,5,11]

    G_2_2[1].names[6]="Vp"  #Vp_0
    G_2_2[1].nodes[6]=3
    G_2_2[1].edges[0:3,6] = [0,8,4]

    G_2_2[1].names[7]="U"  #U_0
    G_2_2[1].nodes[7]=4
    G_2_2[1].edges[0:4,7] = [12,13,9,10]

    # 5/6       19040002000000
    G_2_2[1].weight = 0.250000000000000  
    G_2_2[1].contract_order=[1,4,7,3,5,0,2,6]  

    # OO: [('Vp_0', 1), ('Vp_1', 1), ('V_0', 3), ('V_1', 3)]
    # U_0: [('Up_0', 3), ('oo_2_3', 3), ('V_0', 2), ('V_1', 1)]
    # Up_0: [('Vp_0', 3), ('Vp_1', 2), ('U_0', 1), ('oo_2_3', 1)]
    # V_0: [('Vp_0', 2), ('U_0', 3), ('OO', 3)]
    # V_1: [('U_0', 4), ('oo_2_3', 4), ('OO', 4)]
    # Vp_0: [('OO', 1), ('V_0', 1), ('Up_0', 1)]
    # Vp_1: [('OO', 2), ('Up_0', 2), ('oo_2_3', 2)]
    # oo_2_3: [('Up_0', 4), ('Vp_1', 3), ('U_0', 2), ('V_1', 2)]
    G_2_2[2]=Graphics()
    G_2_2[2].graph_name='G_2_2[2]'
    G_2_2[2].size=8
    G_2_2[2].nodes[:]=-1
    G_2_2[2].edges[:, :]=-1
    G_2_2[2].names[0]="OO"  #OO
    G_2_2[2].nodes[0]=4
    G_2_2[2].edges[0:4,0] = [0,1,2,3]

    G_2_2[2].names[1]="Up"  #Up_0
    G_2_2[2].nodes[1]=4
    G_2_2[2].edges[0:4,1] = [4,5,6,7]

    G_2_2[2].names[2]="V"  #V_0
    G_2_2[2].nodes[2]=3
    G_2_2[2].edges[0:3,2] = [8,9,2]

    G_2_2[2].names[3]="V"  #V_1
    G_2_2[2].nodes[3]=3
    G_2_2[2].edges[0:3,3] = [10,11,3]

    G_2_2[2].names[4]="oo"  #oo_2_3
    G_2_2[2].nodes[4]=4
    G_2_2[2].edges[0:4,4] = [7,12,13,11]

    G_2_2[2].names[5]="Vp"  #Vp_1
    G_2_2[2].nodes[5]=3
    G_2_2[2].edges[0:3,5] = [1,5,12]

    G_2_2[2].names[6]="Vp"  #Vp_0
    G_2_2[2].nodes[6]=3
    G_2_2[2].edges[0:3,6] = [0,8,4]

    G_2_2[2].names[7]="U"  #U_0
    G_2_2[2].nodes[7]=4
    G_2_2[2].edges[0:4,7] = [6,13,9,10]

    # 5/7       6.08040002E+14
    G_2_2[2].weight = 0.250000000000000  
    G_2_2[2].contract_order=[3,4,7,1,5,0,2,6]  

    # OO: [('Vp_1', 1), ('Vp_2', 1), ('V_1', 3), ('V_2', 3)]
    # V_1: [('Vp_1', 2), ('oo_3_4', 3), ('OO', 3)]
    # V_2: [('oo_3_4', 4), ('Vp_2', 3), ('OO', 4)]
    # Vp_1: [('OO', 1), ('V_1', 1), ('oo_3_4', 1)]
    # Vp_2: [('OO', 2), ('oo_3_4', 2), ('V_2', 2)]
    # oo_3_4: [('Vp_1', 3), ('Vp_2', 2), ('V_1', 2), ('V_2', 1)]
    G_2_2[3]=Graphics()
    G_2_2[3].graph_name='G_2_2[3]'
    G_2_2[3].size=6
    G_2_2[3].nodes[:]=-1
    G_2_2[3].edges[:, :]=-1
    G_2_2[3].names[0]="OO"  #OO
    G_2_2[3].nodes[0]=4
    G_2_2[3].edges[0:4,0] = [0,1,2,3]

    G_2_2[3].names[1]="V"  #V_1
    G_2_2[3].nodes[1]=3
    G_2_2[3].edges[0:3,1] = [4,5,2]

    G_2_2[3].names[2]="V"  #V_2
    G_2_2[3].nodes[2]=3
    G_2_2[3].edges[0:3,2] = [6,7,3]

    G_2_2[3].names[3]="Vp"  #Vp_2
    G_2_2[3].nodes[3]=3
    G_2_2[3].edges[0:3,3] = [1,8,7]

    G_2_2[3].names[4]="Vp"  #Vp_1
    G_2_2[3].nodes[4]=3
    G_2_2[3].edges[0:3,4] = [0,4,9]

    G_2_2[3].names[5]="oo"  #oo_3_4
    G_2_2[3].nodes[5]=4
    G_2_2[3].edges[0:4,5] = [9,8,5,6]

    # 4/6       20160004000000
    G_2_2[3].weight = 0.250000000000000  
    G_2_2[3].contract_order=[1,4,0,5,3,2]  


G_3_3 = None
G_3_2 = None
G_22_2 = None
G_22_3 = None



if __name__ == "__main__":
    print G_2_2
    print G_22_2
    print G_22_3
