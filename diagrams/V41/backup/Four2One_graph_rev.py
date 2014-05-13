#!/usr/bin/env python

# First layer |-| |-|

# one 2-site operator

gx1_2_1 = {
    1:[(5,2), (3,3), (8,1)],
    2:[(3,4), (6,3), (8,2)],
    3:[(7,3), (7,4), (1,2), (2,1)],
    4:[(5,3), (6,2), (7,1), (7,2)],
    5:[(8,3), (1,1), (4,1)],
    6:[(8,4), (4,2), (2,2)],
    7:[(4,3), (4,4), (3,1), (3,2)],
    8:[(1,3), (2,3), (5,1), (6,1)]
    }

gx1_2_2 = {
    1:[(5,2), (3,3), (8,1)],
    2:[(3,4), (7,4), (8,2)],
    3:[(4,3), (7,3), (1,2), (2,1)],
    4:[(5,3), (6,2), (3,1), (7,1)],
    5:[(8,3), (1,1), (4,1)],
    6:[(8,4), (4,2), (7,2)],
    7:[(4,4), (6,2), (3,2), (2,2)],
    8:[(1,3), (2,3), (5,1), (6,1)]
}

gx1_2_2 = {
    1:[(5,2), (3,3), (8,1)],
    2:[(3,4), (7,4), (8,2)],
    3:[(4,3), (7,3), (1,2), (2,1)],
    4:[(5,3), (6,2), (3,1), (7,1)],
    5:[(8,3), (1,1), (4,1)],
    6:[(8,4), (4,2), (7,2)],
    7:[(4,4), (6,3), (3,2), (2,2)],
    8:[(1,3), (2,3), (5,1), (6,1)]
}

gx1_2_3 = {
    1:[(3,2), (5,3), (6,1)],
    2:[(5,4), (4,3), (6,2)],
    3:[(6,3), (1,1), (5,1)],
    4:[(6,4), (5,2), (2,2)],
    5:[(3,3), (4,2), (1,2), (2,1)],
    6:[(1,3), (2,3), (3,1), (4,1)]
    }

gx1_2_4 = {
    1:[(7,3), (3,3), (8,1)],
    2:[(3,4), (6,3), (8,2)],
    3:[(7,4), (4,4), (1,2), (2,1)],
    4:[(5,3), (6,2), (7,2), (3,2)],
    5:[(8,3), (7,1), (4,1)],
    6:[(8,4), (4,2), (2,2)],
    7:[(5,2), (4,3), (1,1), (3,1)],
    8:[(1,3), (2,3), (5,1), (6,1)]
}

gx1_3_1 = {
    1:[(5,2), (3,3), (8,1)],
    2:[(3,4), (7,6), (8,2)],
    3:[(7,3), (7,4), (1,2), (2,1)],
    4:[(5,3), (6,1), (7,1), (7,2)],
    5:[(8,3), (1,1), (4,1)],
    6:[(8,4), (4,2), (7,3)],
    7:[(4,3), (4,4), (6,3), (3,1), (3,2), (2,2)],
    8:[(1,3), (2,3), (5,1), (6,1)]
    }

gx1_3_2 = {
    1:[(6,2), (9,4), (10,1)],
    2:[(9,5), (4,3), (10,2)],
    3:[(4,4), (8,3), (10,3)],
    4:[(9,6), (5,4), (2,2), (3,1)],
    5:[(7,3), (8,2), (9,3), (4,2)],
    6:[(10,4), (1,1), (9,1)],
    7:[(10,5), (9,2), (5,1)],
    8:[(10,6), (5,2), (3,2)],
    9:[(6,3), (7,2), (5,3), (1,2), (2,1), (4,1)],
    10:[(1,3), (2,3), (3,3), (6,1), (7,1), (8,1)]
    }

gx1_3_3 = {
    1:[(7,4), (3,3), (8,1)],
    2:[(3,4), (6,3), (8,2)],
    3:[(7,5), (7,6), (1,2), (2,1)],
    4:[(5,3), (6,2), (7,2), (7,3)],
    5:[(8,3), (7,1), (4,1)],
    6:[(8,4), (4,2), (2,2)],
    7:[(5,2), (4,3), (4,4), (1,1), (3,1), (3,2)],
    8:[(1,3), (2,3), (5,1), (6,1)]
    }

gx1_3_4 = {
    1:[(6,2), (4,3), (10,1)],
    2:[(4,4), (9,5), (10,2)],
    3:[(9,6), (8,3), (10,3)],
    4:[(5,4), (9,4), (1,2), (2,1)],
    5:[(6,3), (7,2), (4,1), (9,1)],
    6:[(10,4), (1,1), (5,1)],
    7:[(10,5), (5,2), (9,2)],
    8:[(10,6), (9,3), (3,2)],
    9:[(5,4), (7,3), (8,2), (4,2), (2,2), (3,1)],
    10:[(1,3), (2,3), (3,3), (6,1), (7,1), (8,1)]
    }

gx1_22_1 = {
    1:[(6,2), (4,3), (11,1)],
    2:[(4,4), (10,3), (11,2)],
    3:[(10,4), (8,3), (11,3)],
    4:[(9,3), (9,4), (1,2), (2,1)],
    5:[(6,3), (7,2), (9,1), (9,2)],
    6:[(11,4), (1,1), (5,1)],
    7:[(11,5), (5,2), (10,1)],
    8:[(11,6), (10,2), (3,2)],
    9:[(5,3), (5,4), (4,1), (4,2)],
    10:[(7,3), (8,2), (2,2), (3,1)],
    11:[(1,3), (2,3), (3,3), (6,1), (7,1), (8,1)]
    }

gx1_22_2 = {
    1:[(6,2), (9,3), (11,1)],
    2:[(9,4), (4,3), (11,2)],
    3:[(4,4), (8,3), (11,3)],
    4:[(10,3), (10,4), (2,2), (3,1)],
    5:[(7,3), (8,2), (10,1), (10,2)],
    6:[(11,4), (1,1), (9,1)],
    7:[(11,5), (9,2), (5,1)],
    8:[(11,6), (5,2), (3,2)],
    9:[(6,3), (7,2), (1,2), (2,1)],
    10:[(5,3), (5,4), (4,1), (4,2)],
    11:[(1,3), (2,3), (3,3), (6,1), (7,1), (8,1)]
    }

gx1_22_3 = {
    1:[(6,2), (9,3), (11,1)],
    2:[(9,4), (4,3), (11,2)],
    3:[(4,4), (10,4), (11,3)],
    4:[(5,3), (10,3), (2,2), (3,1)],
    5:[(7,3), (8,2), (4,1), (10,1)],
    6:[(11,4), (1,1), (9,1)],
    7:[(11,5), (9,2), (5,1)],
    8:[(11,6), (5,2), (10,2)],
    9:[(6,3), (7,2), (1,2), (2,1)],
    10:[(5,4), (8,3), (4,2), (3,2)],
    11:[(1,3), (2,3), (3,3), (6,1), (7,1), (8,1)]
    }
    
gx1_22_4 = {
    1:[(9,3), (4,3), (11,1)],
    2:[(4,4), (10,3), (11,2)],
    3:[(10,4), (8,3), (11,3)],
    4:[(9,4), (5,4), (1,2), (2,1)],
    5:[(6,3), (7,2), (9,2), (4,2)],
    6:[(11,4), (9,1), (5,1)],
    7:[(11,5), (5,2), (10,1)],
    8:[(11,6), (10,2), (3,2)],
    9:[(6,2), (5,3), (1,1), (4,1)],
    10:[(7,3), (8,2), (2,2), (3,1)],
    11:[(1,3), (2,3), (3,3), (6,1), (7,1), (8,1)]
    }

gx2_2_1 = {
    1:[(3,3), (3,4), (4,1)],
    2:[(4,2), (3,1), (3,2)],
    3:[(2,2), (2,3), (1,1), (1,2)],
    4:[(1,3), (2,1)]
    }

gx2_2_2 = {
    1:[(3,2), (5,3), (6,1)],
    2:[(5,4), (4,3), (6,2)],
    3:[(6,3), (1,1), (5,1)],
    4:[(6,4), (5,2), (2,2)],
    5:[(3,3), (4,2), (1,2), (2,1)],
    6:[(1,3), (2,3), (3,1), (4,1)]
    }

gx2_3_1 = {
    1:[(5,4), (5,5), (6,1)],
    2:[(5,6), (4,3), (6,2)], 
    3:[(6,3), (5,1), (5,2)],
    4:[(6,4), (5,3), (2,2)],
    5:[(3,2), (3,3), (4,2), (1,1), (1,2), (2,1)],
    6:[(1,3), (2,3), (3,1), (4,1)]
    }

gx2_3_2 = {
    1:[(3,2), (5,4), (6,1)],
    2:[(5,5), (5,6), (6,2)], 
    3:[(6,3), (1,1), (5,1)],
    4:[(6,4), (5,2), (5,3)],
    5:[(3,3), (4,2), (4,3), (1,2), (2,1), (2,2)],
    6:[(1,3), (2,3), (3,1), (4,1)]
    }
