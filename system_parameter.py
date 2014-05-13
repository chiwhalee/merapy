#coding=UTF8
"""
q:
    what is tot_SI_Layer? 

"""
import os

class SystemParam(object):
    param_extra = {}
    @staticmethod
    def precompile_process():
        """ 
        deal with precompile options as in f90 makefile
        """
        SystemParam.MODEL = "Heisenberg"
        #SystemParam.MODEL= "Ising"
        SystemParam.SYMMETRY = "Z2"
        SystemParam.USE_CUSTOM_RANDOM = False
        SystemParam.USE_REFLECTION = False
        SystemParam.NUM_OF_THREADS = None
        SystemParam.only_NN = True
        SystemParam.only_NNN = False

    @classmethod
    def set_default(cls):
        cls.lang_tensor = "py"
        #lang_tensor = "cython"

#mera params        
        #TNet_sys= tensor_network.TensorNetwork()   moved to tensor_network.py
        # below are system default param
        cls.tot_layer = 5   #total layer of mera
        cls.trunc_dim = 2
        cls.nTop = 2
        cls.tot_SI_Layer = 4
        cls.nNeigh = 2
        #      D_0 = 2; D_max=4
#model params
        cls.model_param = {
                "gamma":    1.0, 
                "h":        0.0, 
                "J_NN":     1.0, #see system.cfg there J_NN = -1 
                "J_NNN":    0.10, 
                "alpha":    3.0,  #critical exponent 
                "beta":     -0.10 
                }

#misc
        cls.BackupFile = "UW.dat"
        cls.ReAssumed = False      

    @classmethod
    def system_cfg(cls):
        """
        deal with system.cfg for f90
        """
        #attention_this_may_be_wrong
        #gamma is different from that in system.cfg but...

        cls.model_param_ising= {
                    "gamma":    1.0, 
                    "h":        1.0, 
                    "J_NN":     -1.0, 
                    "J_NNN":    1.0, 
                    "alpha":    2.0,  
                    "beta":     0.0 
                    }
        cls.model_param_heisenberg= {
                    "h":        0.0, 
                    "Jzz":    1.0, 
                    "J_NN":     1.0,  
                    "J_NNN":    0.0, #cls.J_NNN=0.241186
                    "alpha":    2.0,  
                    "beta":     1.0 
                    }

        if cls.MODEL == "Ising": 
            cls.model_param.update(cls.model_param_ising)

            cls.Layer=3
            cls.SI_Layer=2
            cls.D_max=2
        
        elif cls.MODEL == "Heisenberg":
            cls.model_param.update(cls.model_param_heisenberg)
            

            cls.D_max=2
            cls.Layer=3
        else:
            print "error, cls.MODEL is not defined"
            print cls.MODEL

    @classmethod
    def init_all(cls):
        cls.precompile_process()
        cls.set_default()
        cls.system_cfg()
    
    @classmethod
    def show(cls, filename=None):
        import pprint
        #name=["h","gamma","J_NN","J_NNN", "trunc_dim", "tot_layer", "lang_tensor", "MODEL", "USE_REFLECTION"]

        name=["trunc_dim", "tot_layer", "lang_tensor", "SYMMETRY", "MODEL","only_NN", "only_NNN", 
                "NUM_OF_THREADS", "USE_REFLECTION"]
        temp=[(n,cls.__dict__[n]) for n in name]
        #temp.update(cls.param_extra)
        for k in cls.param_extra:
            temp.append((k, cls.param_extra[k]))

        #extra = "OMP_NUM_THREADS"
        #temp += [(extra, os.environ.get(extra))] 
        
        print temp
        #print SystemParam.model_param

        if filename is not None:
            import datetime
            now=str(datetime.datetime.today())
            res= "\n\n!%s\n!%s\n!%s\n"%(now, str(temp), str(SystemParam.model_param))
            out = open(filename, 'a')
            out.write(res)
            out.close()



        #pprint.pprint(temp, indent=4, depth=None )



SystemParam.init_all()



if __name__ == "__main__":
    pass

    #print type(SystemParam.__dict__)
    #print SystemParam.__setattr__("USE_REFLECTION", True)
    #setattr(SystemParam,"USE_REFLECTION",True)
    SystemParam.show()


