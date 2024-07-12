import numpy as np
import time
# from tkinter import *
# from tkinter import ttk
# import matplotlib.pyplot as plt
# from matplotlib.figure import Figure 
# from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
# from tools import render
import subprocess
# import re
HORIZONTAL = 0; VERTICAL = 1

class ADaX():
    """
    \nADaX stands for "AVL designer and executer"
    \nIt allows for automatic generation of .avl files from Python
    \nFunctions preceeded by a '_' take no inputs
    \nFunction preceeded by "__" are private. They are not useful to a user, only to a programmer
    \nSome functions are "code-specific". They help ADaX run, marked by "#code"
    \nSome functions are commands with AVL equivalents, marked by "#avl"
    \nSome functions are identified as custom preset functions, not for general use, marked by "#preset"
    \n┌────────────────────────
    \n│ Text in large boxes
    \n│ represents example code,
    \n│ inputs, or outputs
    \n└────────────────────────
    """
    def __init__(self, planename: str):
        self.__inputlist = ''
        if str.endswith(planename, '.avl'): pass
        else: planename += '.avl'
        self.__planename = planename 
        self.__directory = "AVL/Planes"
        self.__planefile = '{}/{}'.format(self.__directory, planename)
        self.__e1 = [meth for meth in dir(self) if callable(getattr(self,meth))]
        self.__e2 = [getattr(self,meth) for meth in dir(self) if callable(getattr(self,meth))]
        self.__output_list = []

    def output_config(self, output_codes: str|list[str]): #code
        """
        \nExample of an output code:
        \n┌───────────────────────────────────────────
        \n│ 't CLtot' - Total forces,          CLtot
        \n│ 's CYq'   - Stability derivatives, CYq
        \n│ 's CLa'   - Stability derivatives, CLa
        \n│ 'b CXd01' - Body-axis derivatives, CXd01
        \n└───────────────────────────────────────────
        \nPut them into a list to get multiple at once, ex:
        \n┌───────────────────────────────────────────────────────────────    
        \n│ output_codes=['t CLtot', 's CLa']
        \n│ output_codes=['t Alpha', 't CLtot', 't CDind', 't Elevator']
        \n│ output_codes=['s Cma', 's Cnb', 's' Clb']
        \n└───────────────────────────────────────────────────────────────
        \nOutputs will be returned by __output, which is
        \nautomatically called by self.run()
        """
        if type(output_codes) is not list:
            self.__output_list = [output_codes]
        else:
            self.__output_list = output_codes

    def __print(self, text: str): #code
        """
        \nUsed for custom printing of text to Python terminal
        \nSimply used to discern Python text or output text from important messages
        """
        print("[ADaX]", text)

    def __execute(self, execute_code: list[tuple[str]]): #code
        """
        \nGiven a list of execute codes, executes each in order 
        \nExample:
        \n┌────────────────────────────────────────────────────────
        \n│ self.__execute( [('vcv', 'd2 pm 0'), ('vcv', 'a c 0.3')] )
        \n│  →
        \n│ ('vcv', 'd2 pm 0') → self.vcv('d2 pm 0')
        \n│ ('vcv', 'a c 0.3') → self.vcv('a c 0.3')
        \n└────────────────────────────────────────────────────────
        \nUse the exact string of the function name
        \nCan be used to call private functions but why would you do that
        \nCreated for running arbitrary methods in the middle of sweeps, typically vcv()
        """
        print(execute_code)
        for x in execute_code:
            print(x)
            for i,e in enumerate(self.__e1):
                if e == x[0]:
                    self.__e2[i](*x[1:])
                    break
                else: pass

    def __output(self): #code
        """
        \nReturns output_dict in the form
        \n┌─────────────────────────────────
        \n│ {'Alpha': 1.2500
        \n│  'CLtot': 0.2899
        \n│  'CDind': 0.0231 }
        \n└─────────────────────────────────
        \nWhere the keys are determined by 
        \noutput_list (set by output_config())
        \nAutomatically used in self.run()
        """
        if not any(self.__output_list): 
            return {'---':'Output not configured'}
        output_dict = {}
        for out in self.__output_list:
            out = out.split()
            file = self.__get_filename(out[0])
            op = open('AVL/Output/{}'.format(file)).read().split('\n')
            for l in op:
                if f' {out[1]} ' in l:
                    l = l.split()
                    output_dict[out[1]] = l[list.index(l, out[1])+2]
        return output_dict
    
    def __get_filename(self, kword: str): #code
        """
        \nReturns the full string of the output
        \nfile when given a keyword or unique
        \nidentifying string
        \n┌─────────────────────────────────
        \n│ -kword-        -Output file-
        \n│   'b'      Body-axis derivatives
        \n│   's'      Stability derivatives
        \n│   't'         Total forces
        \n└─────────────────────────────────
        \nUsed in __output()
        """
        #Dict for returning proper output file name
        match_key = {'b':'Body-axis derivatives', 's':'Stability derivatives', 't':'Total forces'}
        #If provided kword is in simplest form:
        if kword in 'bst':
            return match_key[kword]
        #If provided kword is truncated form of filename (i don't think this ever happens though?)
        B = 'Body-axis derivatives'
        S = 'Stability derivatives'
        T = 'Total forces'
        logi = [kword in B or kword in B.lower(), kword in S or kword in S.lower(), kword in T or kword in T.lower()]
        if not(logi[0]^logi[1]^logi[2]):
            raise Exception('File code "{}" insufficient to identify unique output file'.format(kword))
        else:
            key = ['b','s','t'][(list.index(logi, True))]
            return match_key[key]

    def input(self, cmd: str): #avl
        """
        \nWrites the "cmd" string to the input list
        \nEquivalent to typing the "cmd" string directly into AVL
        \nThe newline character, "\\n", represents pressing enter
        \nA newline character is automatically added after each call to input()
        """
        # print(str(cmd))
        self.__inputlist += str(cmd) + '\n'

    def vcv(self, vcv: str = 'a a 0'): #avl
        """
        \nVariable, Constraint, Value
        \n
        \nTakes a single string in the same form as the AVL oper menu
        \nCorresponds to driving constraint values in the oper menu
        \n┌───────────────────
        \n│ 'a c 0.42'
        \n│ becomes
        \n│ A lpha -> CL 0.42
        \n└───────────────────
        """
        self.input('{} {} {}'.format(*vcv.split()))

    def _top(self): #avl
        """
        \nTakes you to the top menu
        """
        self.__inputlist += '\n\n\n\n\n\n'

    def _oper(self): #avl
        """
        \nTakes you to the oper menu
        """
        self.__inputlist += '\n\n\n\n\n\noper\n'
        
    def _clear(self): #avl
        """
        \nClears the input list
        """
        self.__inputlist = ''

    def _load(self): #avl
        """
        \nLoads the current plane
        """
        self._top()
        self.input('load {}'.format(self.__planefile))

    def _save(self) -> None: #avl
        """
        \nSaves all three outputs
        """
        self.input('st AVL/Output/Stability derivatives')
        self.input('o')
        self.input('ft AVL/Output/Total forces')
        self.input('o')
        self.input('sb AVL/Output/Body-axis derivatives')
        self.input('o')
        self.input('fn AVL/Output/Surface forces')
        self.input('o')
        self.input('fs AVL/Output/Strip forces')
        self.input('o')


    def _x(self) -> None: #avl
        """
        \nExecute run case you set up 
        """
        self.input('x')
    
    def use_run(self, run_filename:str, case_number:int) -> None:
        """
        \nLoads an existing .run file and chooses a case
        \nTakes the .run file name as it appears in "AVL/Planes"
        \nTakes case_number as an integer
        \nReturns nothing
        """
        run_filename = run_filename.strip()
        if not run_filename.endswith('.run'): run_filename += '.run'
        self.input('f')
        self.input('AVL/Planes/{0}'.format(run_filename))
        self.input(case_number)
        self.input('c1\n')

    def run(self, print_output:bool=True) -> dict[str:str]: #avl
        """
        \nRun AVL and input all commands from input list
        """
        self._top()
        self.input('quit')
        self.sp = subprocess.Popen('avl.exe',
                                   shell=False,
                                   stdin=subprocess.PIPE,
                                   stdout=open('Log.txt', 'w'), 
                                   stderr=subprocess.PIPE)
        self.sp.stdin.write(self.__inputlist.encode('utf-8'))
        self.sp.stdin.flush()
        self.sp.communicate()
        self._clear()
        out = self.__output()
        if print_output: 
            self.__print("Run results:")
            for key, item in out.items():
                print('{0:7}{1:<9} {2:<12}'.format(' ',key,item))
        return out
    
    def sweep(self, vc: str, v_array: list, executes: str|list[str]): #preset
        """
        \n┌────────────────────────────────────────────────
        \n│  -input-             -description-
        \n│    vc       Variable-constraint pair for vcv()
        \n│  v_array         Value array for vcv()
        \n│  executes          Additional executes
        \n└────────────────────────────────────────────────
        \nReturns a dictionary of output dictionaries 
        \n'outs' in the form of
        \n┌────────────────────────────────────────────────
        \n│ outs[vc + v_array[ 0]] = output_dict  0
        \n│ outs[vc + v_array[ 1]] = output_dict  1
        \n│         ...
        \n│ outs[vc + v_array[-1]] = output_dict -1
        \n└────────────────────────────────────────────────
        \nRoom for improvement: Instead of opening and running AVL
        \nfor each vcv, just run then all in one instance of AVL.
        \nYou'll need a buffer folder holding every vcv's results,
        \nbut you might get better speed at the cost of a folder with
        \n100+ files in it (delete the files after reading them ofc)
        """
        self.__print(("Beginning sweep over {} {} to {} {}".format(vc, v_array[0], vc, v_array[-1])))
        adsub = ADaX(self.__planename)
        adsub.output_config(self.__output_list)
        outs = {}
        header = ['case no.', 'sweep vcv', 'case time (s)', 'est. time remaining (s)']
        print('{:>10}  {:>10}  {:>15}  {:>20}'.format(*header))
        print('─'*70)
        n = len(v_array)
        t0 = time.time()
        dt = np.array([], dtype=int)
        for i,v in enumerate(v_array):
            t1 = time.time()
            adsub._load()
            adsub._oper()
            adsub.__execute(executes)
            adsub.vcv(vc + ' {}'.format(v))
            adsub._x()
            adsub._save()
            outs[vc + ' {}'.format(v)] = adsub.run(print_output=False)
            t2 = time.time()
            dt = np.append(dt, t2 - t1)
            print('{:>10}  {:>10}  {:>15}  {:>20}'.format(f'{i+1}/{n}', vc+f' {v}', np.around(t2-t1,2), np.around(self.__time_remaining(dt,i,n),2)))
        self.__print("Finished sweep with total time {:<12.2f}".format(time.time()-t0))
        return outs
    
    def __time_remaining(self, dt: np.ndarray[int], i: int, n: int) -> int: #code
        """
        \ndt: case time array
        \ni: current case
        \nn: total cases
        \nUsed in sweep() to show the remaining time
        """
        i += 1
        n_remaining = n-i
        if n_remaining == 0:
            return 0.0
        if len(dt) == 1:
            return dt[0]*n_remaining
        ia = np.arange(1,float(i)+0.1,1)
        iz = np.ones_like(ia)
        a = np.stack([ia, iz], axis=1)
        c = np.linalg.lstsq(a,dt,rcond=None)[0]
        return n_remaining*c[1] + i*c[0]      
        
    def inviscid_polar(self, a_array: np.ndarray, trim: str) -> None: #preset
        """
        \nGenerates a drag polar file for a range of alpha
        \n⚠ Name your elevator "Elevator" for best results
        \nTakes a_array as a numpy array of desired angles of attack
        \nTakes a trim vcv, should be 'd3 pm 0' or something, where d3 is the elevator
        \nDoes not return anything, generates "AVL/Output/Invis Polar.dat"
        \n
        \nRoom for improvement: write results to "Invis Polar.dat" after each run;
        \ncurrently it attempts to write it all at the end. If there's any error,
        \nnothing gets written.
        """
        self.output_config(['t Alpha', 't CLtot', 't CDind', 't Elevator'])
        polar_dict = self.sweep(vc='a a', v_array=a_array, executes=[('vcv', trim)])
        f = open('AVL/Output/Invis Polar.dat', 'w')
        f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format('Alpha', 'CLtot', 'CDind', 'Elevator'))
        f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format('-'*20, '-'*20, '-'*20, '-'*20))
        for k, val in polar_dict.items():
            f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format(val['Alpha'], val['CLtot'], val['CDind'], val['Elevator']))

    def stall_prediction(self, Clx: float) -> float: #preset
        """
        \nStall prediction based on airfoil maximum Cl (Clx)
        \nChecks strip forces along wing
        \nMakes the assumption that if strip Cl reaches 95% of Clx, wing stalls
        \n⚠ AVL is a potential flow solver; it cannot accurately predict stall! ⚠
        """
        Clx = 0.95*Clx
        adsub = ADaX(planename=self.__planename)
        adsub.output_config(self.__output_list)
        a = 0
        k = -10
        Cl_error = -1
        self.__print("Beginning stall prediction iteration with Clx = {}".format(Clx))
        header = ['iter no.', 'Alpha', 'Max strip Cl', 'res']
        print('{:>10}  {:>10}  {:>13}  {:>10}'.format(*header))
        print('─'*55)
        n = 0
        while abs(Cl_error) > 1e-3:
            a += Cl_error*k
            n += 1
            adsub._load()
            adsub._oper()
            adsub.vcv('a a {}'.format(a))
            adsub._x()
            adsub._save()
            adsub.run(print_output=False)
            content = open('AVL/Output/Strip forces').readlines()
            cl = []
            for i, line in enumerate(content):
                if ' j ' in line:
                    for j,line in enumerate(content[i+1:]):
                        if line == '\n': 
                            break
                        else: 
                            cl.append(line.split()[9])
                    cl = max(np.array(list(map(float,cl))))
                    break
            Cl_error = cl - Clx
            header = [n, a, cl, Cl_error]
            print('{:>10}  {:>10.4f}  {:>13.4f}  {:>10.5f}'.format(*header))
        self.__print("Finished stall prediction at Alpha = {}".format(np.around(a,4)))

class des_obj():
    """
    \nParent class for design objects, since they all have a similar structure.
    \nContains functions common to all design objcets; does nothing on its own.
    """

    def __init__(self):
        pass

    def modify(self, **kwargs):
        "Modify an existing value in .config"
        keys = self.config.keys()
        for kw, val in kwargs.items():
            if kw in keys: 
                self.config[kw] = val
            else: 
                raise KeyError('Cannot modify {}: KeyError "{}"'.format(self, kw))
    
    
    def assign(self, **kwargs):
        "Only use this to append non-standard parameters to des objects"
        for kw, val in kwargs.items():
            self.config[kw] = val
        self.format += '1'
        self.kw += '1'

    
    def __assign(self, kwargs):
        "Used for __init__ assignments"
        keys = self.config.keys()
        for kw, val in kwargs.items():
            if kw in keys: 
                self.config[kw] = val
            else: 
                raise KeyError('Cannot modify {}: KeyError "{}"'.format(self, kw))

    
    def __lv_pair(self, keys: str | list[str], cmt: bool = True):
        "Print label-value pairs to .avl files"
        top = ''
        bot = ''
        if type(keys) is str:
            keys = [keys]
        for i in keys:
            if cmt:
                top += '{:<12}  '.format('#{}'.format(i))
            else:
                top += '{:<12}  '.format('{}'.format(i)) 
            val = self.config[i]
            if type(val) is str: bot += '{:<12}  '.format(val)
            elif type(val) is list or type(val) is tuple: 
                for v in val:
                    v = np.around(v,5)
                    bot += '{:<3} '.format(v)
                    if i == 'TRANSLATE': bot += ' '
                bot += '  '
            elif val is None:
                return ''
            else: 
                val = np.around(val,5)
                bot += '{:<12}  '.format(val)
        return top + '\n' + bot + '\n'
    
    
    def config_to_avl(self, f):
        "Handles des object .config print order"
        keys = []
        for i, label in enumerate(self.config.keys()):
            keys.append(label)
            if self.format[i] == '1':
                if self.kw[i] == '1': cmt = False
                else: cmt = True
                f.write(self.__lv_pair(keys=keys, cmt=cmt))
                keys = []
            else: continue

    
    def children(self, *args):
        "Print all children of a des object"
        print('Children of {}:'.format(self))
        try: children = self.surfs
        except: pass
        try: children = self.secs
        except: pass
        try: children = self.ctrls
        except: pass
        if type(children) is list:
            for i, ch in enumerate(children):
                print(' {}. {}'.format(i, ch))
        else:
            print(' 1. {}'.format(children))

    
    def parents(self, *args):
        "Print all parents of a des object"
        print('Parent of {}:'.format(self))
        print('    {}'.format(self.parent))
    
class des_plane(des_obj):
    """
    \nPlane/header component of .avl files
    \nParent: None
    \nChildren: des_surf
    """

    def __init__(self, **kwargs):
        #AVL equivalent parameters
        self.config = {'Name': 'Unnamed plane',
                       'Mach': 0.0,
                       'iYsym': (0), 'iZsym': int(0), 'Zsym': 0.0,
                       'Sref': 0.0, 'cref': 0.0, 'bref': 0.0,
                       'Xref': 0.0, 'Yref': 0.0, 'Zref': 0.0,
                       'CDp': 0.0}
        #Formatting code for .avl file
        self.format = '110010010011'
        self.kw = '000000000000'
        #Keep track of child surfaces
        self.surfs = []
        #Config assignment
        self._des_obj__assign(kwargs)
        #If custom surface is defined
        self.__reference = None
        #Rendering variables
        self.surf_book = {}
        self.ctrl_book = {}

    
    def __str__(self):
        """
        Custom printing
        """
        return 'des_plane object "{}"'.format(self.config['Name'])
    
    
    def write_to_avl(self, filename: str):
        """
        Begin writing the geometry to .avl file
        """
        print('[Geometry] Begin writing geometry to AVL/Planes/{0}.avl'.format(filename))
        width = 95
        if self.__reference:
            print('[Geometry] Using reference surface {} to update reference values'.format(self.__reference))
            self.__reference.calculate_reference_values()
            self.assign(Sref=self.__reference.Sref, bref=self.__reference.bref, cref=self.__reference.Cref)
        f = open('AVL/Planes/{0}.avl'.format(filename), 'w')
        self.config_to_avl(f)
        for sf in self.surfs:
            f.write('#' + '='*width + '\n#' + '='*width + '\nSURFACE\n')
            sf.config_to_avl(f)
            for sc in sf.secs:
                f.write('#' + '-'*width + '\nSECTION\n')
                sc.config_to_avl(f)
                for ct in sc.ctrls:
                    f.write('\nCONTROL\n')
                    ct.config_to_avl(f)
        print('[Geometry] Finished writing geometry to AVL/Planes/{0}.avl'.format(filename))

    def use_reference(self, reference_surface):
        """
        \nUses a des_surf as reference surface
        \nReference quantites are calculated when nwrite_to_avl() is called
        """
        self.__reference=reference_surface

class des_surf(des_obj):
    """
    \nSurface component of .avl files
    \nParent: des_plane
    \nChildren: des_sec
    """

    def __init__(self, parent: des_plane, **kwargs):
        #AVL equivalent parameters
        self.config = {'Name': 'Unnamed surface',
                       'Nchordwise': 8, 'Cspace': 1.0, 'Nspanwise': 12, 'Sspace': 1.0,
                       'YDUPLICATE': None,
                       'TRANSLATE': (0.0,0.0,0.0),
                       'CDCL': None}
        #Formatting code for .avl file
        self.format = '10001111'
        self.kw = '00000111'
        #Parent-child coupling
        self.parent = parent
        self.parent.surfs.append(self)
        self.secs = []
        #Assign the stuff
        self._des_obj__assign(kwargs)
        #Internal variables
        self.reference = []
        self.Sref = 0
        self.Cref = 0
        self.bref = 0

    
    def __str__(self):
        """
        Custom printing
        """
        return 'des_surf object "{}"'.format(self.config['Name'])

    
    def calculate_reference_values(self):
        """
        Calculate surface reference values - only used for custom reference
        """
        Yle = [sc.config['Yle'] for sc in self.secs]
        Chord = np.array([sc.config['Chord'] for sc in self.secs])
        c = np.array([])
        y = np.array([])
        for i in range(len(Yle)-1):
            c = np.append(c, np.linspace(Chord[i], Chord[i+1], 101))
            y = np.append(y, np.linspace(Yle[i], Yle[i+1], 101))
        self.Sref = 2*np.trapz(c,y)
        self.Cref = 2*np.trapz(c**2, y)/self.Sref
        self.bref = 2*Yle[-1]

    
    def order(self, orientation: int = 0):
        """
        \nSort surface sections into proper order
        \n0 horizontal
        \n1 vertical
        """
        coord_of_interest = ['Yle','Zle'][orientation]
        Yle_collection = []
        for i, sc in enumerate(self.secs):
            Yle_collection.append([sc.config[coord_of_interest], i, sc])
        Yle_sorted = sorted(Yle_collection, key=lambda x: x[0])
        for i, new_sc in enumerate(Yle_sorted):
            self.secs[i] = new_sc[2]

class des_sec(des_obj):
    """
    \nSection component of .avl files
    \nParent: des_surf
    \nChildren: des_ctrl
    """

    def __init__(self, parent: des_surf, **kwargs):
        #AVL equivalent parameters
        self.config = {'Xle': 0.0, 'Yle': 0.0, 'Zle': 0.0, 'Chord': 1.0,
                       'Ainc': 0.0, 
                       'Nspanwise': 0, 'Sspace': 0,
                       'AFILE': None}
        #Formatting code for .avl file
        self.format = '00000011'
        self.kw = '00000001'
        #Parent-child coupling
        self.parent = parent
        self.parent.secs.append(self)
        self.ctrls = []
        #Assign that stuff
        self._des_obj__assign(kwargs)

    
    def __str__(self):
        """
        Custom printing
        """
        TRANSLATE = self.parent.config['TRANSLATE']
        if TRANSLATE is None:
            TRANSLATE = (0,0,0)
        else: pass
        global_pos = (self.config['Xle'] + TRANSLATE[0], 
                      self.config['Yle'] + TRANSLATE[1], 
                      self.config['Zle'] + TRANSLATE[2])
        return 'des_sec object located at (Xle, Yle, Zle) = ({:> 7g} {:> 7g} {:> 7g})'.format(*np.around(global_pos,3))
    
    
    def interp(s1, s2, span: float = 0, orientation: int = 0):
        """
        Create a new section using linear interpolation at a span between two existing sections
        0 horizontal
        1 vertical
        """
        spanwise_variable = ['Yle', 'Zle'][orientation]
        not_the_spanwise_variable = ['Zle', 'Yle'][orientation]
        val = []
        x1 = [s1.config[spanwise_variable], s2.config[spanwise_variable]]
        for kw in ['Xle', not_the_spanwise_variable, 'Chord', 'Ainc']:
            x2 = [s1.config[kw], s2.config[kw]]
            val.append(np.interp(span, x1, x2))
        if orientation:
            return des_sec(parent=s1.parent, Xle=val[0], Yle=val[1], Zle=span, Chord=val[2], Ainc=val[3], AFILE=s1.config['AFILE'])
        else: 
            return des_sec(parent=s1.parent, Xle=val[0], Yle=span, Zle=val[1], Chord=val[2], Ainc=val[3], AFILE=s1.config['AFILE'])

class des_ctrl(des_obj):
    """
    \nControl component of .avl files
    \nParent: des_sec
    \nChildren: None
    """

    def __init__(self, parent: des_sec,  **kwargs):
        #AVL equivalent parameters
        self.config = {'Cname': 'Unnamed control', 'Cgain': 1.0, 'Xhinge': 0.7, 'HingeVec': (0.0, 0.0, 0.0), 'SgnDup': 1.0}
        #Formatting code for .avl file
        self.format = '00001'
        self.kw = '00000'
        #Parent-child coupling
        self.parent = parent
        self.parent.ctrls.append(self)
        #Assign the stuff
        self._des_obj__assign(kwargs)

    def __str__(self):
        """
        Custom printing
        """
        return 'des_ctrl object "{}"'.format(self.config['Cname'])

class Case_Maker:
    """
    \nClass for handling custom case creation
    \nInit -> open existing -> insert new run -> write to .run
    """
    def __init__(self, avl_file:str, run_file:str):
        """
        \nAssigns plane file and run file locations
        """
        avl_file = avl_file.strip()
        run_file = run_file.strip()
        if not avl_file.endswith('.avl'): avl_file += '.avl'
        if not run_file.endswith('.run'): run_file += '.run'
        self.__planename = avl_file
        self.__runfile = run_file

    def _open_existing_dot_run(self) -> None:
        """
        \nOpens the existing run file and reads existing run cases
        \nStores in self.__existing_run. Keys are case names
        \nIf no existing file, self.__existing_run is empty.
        """
        self.__existing_run = {}
        try:
            contents = open('AVL/Planes/{}'.format(self.__runfile), 'r').readlines()
        except: 
            return
        line_break_indices = []
        for i,line in enumerate(contents):
            if ' ---------------------------------------------\n' in line:
                line_break_indices.append(i)
        line_break_indices.append(i+2)
        for break_no,break_index in enumerate(line_break_indices[:-1]):
            #Stephen King wishes he could match how horrifying this line is
            self.__existing_run[contents[break_index+1].split(':  ')[1].strip()] = contents[break_index-1:line_break_indices[break_no+1]-1] 
        # print(self.__existing_run)

    def insert_new_run_case(self, case_name:str, vcv:list[str], flight:dict[str:float], mass:dict[str:float], viscous:dict[str:float]) -> int:
        """
        \nInserts new run case into self.__existing_run
        \n- Prompts user to overwrite if case name already exists
        \n- Retrieves default values (calls self.fetch_defaults())
        \n- Converts vcv to run case format (calls self.vcv_to_run_case())
        \n- If overwriting, replaces existing run #. If not, appends to end of run
        \nDoes NOT write to .run file
        \nReturns 100 if user approves overwriting or insert is successful
        \nReturns 200 if user denies overwriting
        """
        #Determine the case number
        case_number = len(self.__existing_run) + 1
        if case_name in self.__existing_run.keys():
            print("[CaseMaker] Case {0} already exists in {1}. Overwrite? 1/0 = Y/N".format(case_name, self.__runfile))
            while True:
                user_input = input("Enter:")
                if user_input in ('1','0'):
                    break
            if user_input: 
                case_number = float(self.__existing_run[case_name][2].split(':  ')[0].split()[-1])
            else: return 200
        #Retrieve any default values (marked as -1)
        igor_fetch_me_the_defaults = []
        for key,item in flight.items():
            if item == -1: igor_fetch_me_the_defaults.append(key)
        for key,item in mass.items():
            if item == -1: igor_fetch_me_the_defaults.append(key)
        if any(igor_fetch_me_the_defaults):
            basket = self.fetch_defaults(shopping_list=igor_fetch_me_the_defaults)
            for key,item in flight.items():
                if key in igor_fetch_me_the_defaults: flight[key] = basket.pop(0)
            for key,item in mass.items():
                if key in igor_fetch_me_the_defaults: mass[key] = basket.pop(0)
        #Convert symbol vcv to text vcv (e.g., 'a' -> "Alpha")
        converted_vcv = self.vcv_to_run_case(vcv_list=vcv)
        #Create the contents of the new run case as a list
        new_case_content = [
            "\n"
        ,   " ---------------------------------------------\n"
        ,   " Run case {0:>2g}:  {1:<}\n".format(case_number, case_name)
        ,   "\n"
        ]
        for single_vcv in converted_vcv:
            var, con, val = single_vcv.split()
            val = float(val)
            new_case_content.append(" {0:<13}->  {1:<12}=  {2:< 12.5E}\n".format(var, con, val))
        #Format is not automatically adjusted for value magnitude. 6-sig. fig.-Scientific notation is used for all.
        default_format = " {0:<10}=  {1:< 12.5E} {2}\n"
        units = self.__units_config()
        flight_parameters = []
        mass_parameters = []
        viscous_parameters = []
        #I guess you could replace this with a function
        for key, item in flight.items():
            try: flight_parameters.append(default_format.format(key, item, units[key]))
            except: flight_parameters.append(default_format.format(key, item, ''))
        flight_parameters.append("\n")
        for key, item in mass.items():
            try: mass_parameters.append(default_format.format(key, item, units[key]))
            except: mass_parameters.append(default_format.format(key, item, ''))
        for key, item in mass.items():
            try: viscous_parameters.append(default_format.format(key, item, units[key]))
            except: viscous_parameters.append(default_format.format(key, item, ''))
        new_case_content.extend(flight_parameters)
        new_case_content.extend(mass_parameters)
        new_case_content.extend(viscous_parameters)
        self.__existing_run[case_name] = new_case_content
        return 100
    
    def _write_the_dot_run_file(self) -> None:
        """
        \nWrites the contents of self.__existing_run to a .run file
        \nTakes nothing, returns nothing
        """
        print("[CaseMaker] Begin writing AVL/Planes/{}".format(self.__runfile))
        with open('AVL/Planes/{}'.format(self.__runfile), 'w') as f:
            for existing_case_name, existing_case_list in self.__existing_run.items():
                for line in existing_case_list:
                    f.write(line)
        print("[CaseMaker] Finished writing AVL/Planes/{}".format(self.__runfile))
        
        
            
    def vcv_to_run_case(self, vcv_list:list[str]) -> list[str]:
        """
        \nConverts the vcv to the word version used by the .run files
        \nTakes a list[str] of vcv
        \nReturns a converted list[str] where the vcv symbols are replaced with the word versions
        """
        alias = {
            'a' : 'alpha'
        ,   'b' : 'beta'
        ,   'r' : 'pb/2V'
        ,   'p' : 'qc/2V'
        ,   'y' : 'rb/2V'
        ,   'c' : 'CL'
        ,   's' : 'CY'
        ,   'rm': 'Cl roll mom'
        ,   'pm': 'Cm pitchmom'
        ,   'ym': 'Cn yaw  mom'
        }
        self._fetch_control_names()
        for d, name in self.__control_names.items():
            alias[d] = name
        converted_list = []
        for vcv in vcv_list:
            var, con, val = vcv.split()
            try: var = alias[var]
            except: print("[CaseMaker] ⚠ Variable {} not aliased.".format(var))
            try: con = alias[con]
            except: print("[CaseMaker] ⚠ Constraint {} not aliased.".format(con))
            converted_list.append("{} {} {}".format(var,con,val))
        return converted_list
    
    def fetch_defaults(self, shopping_list: list[str]) -> list[str]:
        """
        \nGets default case values from AVL
        \nShopping list is a list of strings with the names of the desired parameters
        \ne.g.: shopping_list = ['X_cg', 'Y_cg', 'Z_cg', 'mass', 'alpha', 'beta']
        \nReturns values in the same order they are requested
        """
        print('[CaseMaker] Fetching case defaults from AVL')
        ad = ADaX(planename=self.__planename)
        ad._load()
        ad._oper()
        ad._x()
        ad.input('s')
        ad.input('AVL/Output/Default run.run')
        ad.input('y')
        ad.run(print_output=False)
        contents = open('AVL/Output/Default run.run').readlines()
        basket = []
        for item in shopping_list:
            for line in contents:
                if item in line:
                    suspect = line.strip().split('=')[1]
                    if ' ' in suspect: basket.append(float(suspect.split()[0]))
                    else: basket.append(float(suspect))
        print('[CaseMaker] Fetch complete')
        return basket
        
    def _fetch_control_names(self) -> dict[str:str]:
        """
        \nGets your control surface names from AVL
        \nCreates self.__control_names, a dictionary of form dict[str:str] 
        \ne.g.: {"d1":"Flaps", "d2":"Aileron", "d3":"Elevator", "d4":"Rudder"}
        \nReturns the dictionary for potential future development
        """
        print('[CaseMaker] Fetching control names from AVL')
        ad = ADaX(planename=self.__planename)
        ad._load()
        ad._oper()
        ad._x()
        ad._save()
        ad.run(print_output=False)
        contents = open('AVL/Output/Total forces').read()
        contents = contents.split('Trefftz')[1].split('\n')[3:]
        i = 0
        line = contents[i]
        self.__control_names = {}
        while line != '':
            i += 1
            self.__control_names['d{0}'.format(i)] = line.split('=')[0].strip()
            line = contents[i]
            # print(self.__control_names)
        print('[CaseMaker] Fetch complete')
        return self.__control_names

    def __units_config(self) -> dict[str:str]:
        """
        \nReturns a dict[str:str] of units for 
        """
        units = {
            "alpha":"deg"
        ,   "beta":"deg"
        ,   "bank":"deg"
        ,   "elevation":"deg"
        ,   "heading":"deg"
        ,   "velocity":"Lunit/Tunit"
        ,   "density":"Munit/Lunit^3"
        ,   "grav.acc":"Lunit/Tunit^2"
        ,   "turn_rad":"Lunit"
        ,   "X_cg":"Lunit"
        ,   "Y_cg":"Lunit"
        ,   "Z_cg":"Lunit"
        ,   "mass":"Munit"
        ,   "Ixx":"Munit-Lunit^2"
        ,   "Iyy":"Munit-Lunit^2"
        ,   "Izz":"Munit-Lunit^2"
        ,   "Ixy":"Munit-Lunit^2"
        ,   "Iyz":"Munit-Lunit^2"
        ,   "Izx":"Munit-Lunit^2"
        }
        return units
        
class tools:
    """
    \nExtra tools used for doing fun stuff I really don't know how else to explain it
    """

    def __init__(self):
        pass

    def render(pts: np.ndarray, beta: float, gamma: float) -> np.ndarray:
        """
        \nRotate points from aircraft body reference to render in 3D
        \nGive pts in aircraft body axis and beta/gamma in degrees
        \nReturns points as (X,Y) viewed from beta/gamma 
        \nRunning plt.plot(X,Y) will display the rendered points
        \n
        \n---Example---
        \npt = np.array([[0.0,1.0,1.1,0.4,0.0,1.0,1.1,0.4,0.0], [0.0,0.0,2.4,2.4,0.0,0.0,-2.4,-2.4,0.0], [0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0]])
        \nb = 25
        \ng = 45
        \nrotated = render(pts = pt, gamma=g, beta=b)
        \nplt.plot(rotated[0], rotated[1])
        \nplt.axis('equal')
        \nplt.show()
        \n-------------
        """
        b = beta*np.pi/180
        sb = np.sin(b)
        cb = np.cos(b)
        g = gamma*np.pi/180
        sg = np.sin(g)
        cg = np.cos(g)
        C2 = np.array([[cg*cb, -sg, cg*sb],[sg*cb, cg, sg*sb],[-sb, 0, cb]])
        Xhat = np.array([[0],[-1],[0]])
        Yhat = np.array([[0],[0],[1]])
        Xhat = np.matmul(C2,Xhat)
        Yhat = np.matmul(C2,Yhat)
        P = np.array([*np.transpose(Xhat),*np.transpose(Yhat)]) 
        return np.matmul(P,pts)
    
    def _autodebug():
        """
        \nLooks for error messages in "Log.txt" to relay problems to user
        """
        print("[Autodebug] Begin debug scan of Log.txt")
        error_book =  {"corrupted": "Suspected error finding or loading .avl file. Double check directory or file name.",
                       "first!": "Failed to access results; flow execution not completed. Is _x() present in __AVL_Commands__?",
                       "failed": "Trim convergence failed. Double check which control surface is trimming.",
                       "large": "Resulting alpha is too large.",
                       "zero-camber": "Airfoil incorrectly loaded. Check the airfoil name, directory, or existence.",
                       "not recognized": "Unrecognized command passed to AVL. Check spelling and which menu you are in.",
                       "    **": "Conflicting vcv.",
                       "****": "NaN result. If you've managed this, please take a screenshot cause I want to know how.",
                       "Insufficient": "Insufficient spanwise points. Try increasing SURFACE NSpanwise or defining SECTION NSpanwise"}
        log = open("Log.txt").readlines()
        for flag_term in error_book.keys():
            for i, line in enumerate(log):
                if flag_term in line:
                    print("[Autodebug]", "Line {0}:".format(i), error_book[flag_term])
                else: pass
        print("[Autodebug] Finished debug scan")
        