
import argparse, argcomplete

if 1:     
    parser = argparse.ArgumentParser(description='parser for meassure and analysis')
    parser.add_argument('dir_list', nargs='*', default=['./'])
    #parser.add_argument('-w', '--which', nargs='*', default=None)
    parser.add_argument('-ew', '--exclude_which', nargs='*', default=None)
    parser.add_argument('-d', '--dim_list', nargs='*',  type=int,  default=None)
    parser.add_argument('-l', '--layer_list', nargs='*', type=int, default=None)
    parser.add_argument('-sh', '--mera_shape_list', nargs='*', default='all')
    parser.add_argument('-sh_min', '--sh_min', nargs=2, type=int,  default=None)
    parser.add_argument('-sh_max', '--sh_max', nargs=2, type=int, default=None)
    parser.add_argument('-f', '--force', action='store_true', default=False) 
    parser.add_argument('-ft', '--fault_tolerant', type=int, default=1) 
    parser.add_argument('-r', '--recursive', action='store_true', default=False) 
    
   
    #parser.add_argument('-f', '--force', action='store_true')
    if 0: 
        args = parser.parse_args()
        args = vars(args)
    if 0: 
        if not args['mera_shape_list'] in [None, 'all']: 
            dim_list = args['dim_list']; args.pop('dim_list')
            layer_list = args['layer_list']; args.pop('layer_list')
            mera_shape_list = [(d, l) for d in dim_list for j in layer_list]
            args['mera_shape_list'] = mera_shape_list
        
        
if __name__ == '__main__': 
    
    args = parser.parse_args()
    args = vars(args)
    print args

