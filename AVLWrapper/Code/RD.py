import numpy as np
import time
# from tkinter import *
# from tkinter import ttk
import matplotlib.pyplot as plt
# from matplotlib.figure import Figure 
# from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk) 
# from tools import render
import subprocess
import re
HORIZONTAL = 0; VERTICAL = 1

class RAMDaX():
    """
    \n══════════
    \nUSER MANUAL
    \n══════════
    \nRAMDaX stands for ".run, .avl,. mass designer and executer"
    \nIt allows for automatic and manual generation of .run, .avl, and .mass files from Python
    \nFunctions preceeded by a '_' take no inputs
    \nFunction preceeded by "__" are private. They are not useful to a user, only to a programmer
    \nSome functions are "code-specific". They help RAMDaX run, marked by "#code"
    \nSome functions are commands with AVL equivalents, marked by "#avl"
    \nSome functions are identified as custom preset functions, not for general use, marked by "#pre"
    \n┌────────────────────────
    \n│ Text in large boxes
    \n│ represents example code,
    \n│ inputs, or outputs
    \n└────────────────────────
    """
    def __init__(self, planename: str):
        self.__inputlist = ''
        self.__directory = "AVL/Planes"
        if str.endswith(planename, '.avl'): pass
        else: planename += '.avl'
        self.__planename = planename 
        self.__planefile = '{}/{}'.format(self.__directory, planename)
        self.e1 = [meth for meth in dir(self) if callable(getattr(self,meth))]
        self.e2 = [getattr(self,meth) for meth in dir(self) if callable(getattr(self,meth))]
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
        \nautomatically called by _run
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
        print("[RD]", text)

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
        for x in execute_code:
            for i,e in enumerate(self.e1):
                if e == x[0]:
                    self.e2[i](*x[1:])
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
        \nAutomatically used in _run()
        """
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
        def match_key(kword): #c
            match kword.lower():
                case 'b':
                    return 'Body-axis derivatives'
                case 's':
                    return 'Stability derivatives'
                case 't':
                    return 'Total forces'
                case _:
                    pass
        if kword in 'bst':
            return match_key(kword)
        B = 'Body-axis derivatives'
        S = 'Stability derivatives'
        T = 'Total forces'
        logi = [kword in B or kword in B.lower(), kword in S or kword in S.lower(), kword in T or kword in T.lower()]
        if not(logi[0]^logi[1]^logi[2]):
            raise Exception('File code "{}" insufficient to identify unique output file'.format(kword))
        else:
            key = ['b','s','t'][(list.index(logi, True))]
            return match_key(key)

    def input(self, cmd: str): #avl
        """
        \nWrites the "cmd" string to the input list
        \nEquivalent to typing the "cmd" string directly into AVL
        \nThe newline character, "\\n", represents pressing enter
        \nA newline character is automatically added after each call to input()
        """
        self.__inputlist += cmd + '\n'

    def vcv(self, vcv: str = 'a a 0'): #avl
        """
        \nVariable, Constraint, Value
        \n
        \nTakes a single string in the same form as the AVL oper menu
        \nCorresponds to driving constraint values in the oper menu
        \n┌───────────────────
        \n│ 'a c 0.42'
        \n│ →
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

    def _save(self): #avl
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


    def _x(self): #avl
        """
        \nExecute run case you set up 
        """
        self.input('x')

    def _run(self): #avl
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
        # print(out)
        return out
    
    def sweep(self, vc: str, v_array: list, executes: str|list[str]): #pre
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
        """
        self.__print(("Beginning sweep over {} {} to {} {}".format(vc, v_array[0], vc, v_array[-1])))
        rdsub = RAMDaX(self.__planename)
        rdsub.output_config(self.__output_list)
        outs = {}
        header = ['case no.', 'sweep vcv', 'case time (s)', 'est. time remaining (s)']
        print('{:>10}  {:>10}  {:>15}  {:>20}'.format(*header))
        print('─'*70)
        n = len(v_array)
        t0 = time.time()
        dt = np.array([], dtype=int)
        for i,v in enumerate(v_array):
            t1 = time.time()
            rdsub._load()
            rdsub._oper()
            for x in executes:
                rdsub.__execute(execute_code=x)
            rdsub.vcv(vc + ' {}'.format(v))
            rdsub._x()
            rdsub._save()
            outs[vc + ' {}'.format(v)] = rdsub._run()
            t2 = time.time()
            dt = np.append(dt, t2 - t1)
            print('{:>10}  {:>10}  {:>15}  {:>20}'.format(f'{i+1}/{n}', vc+f' {v}', np.around(t2-t1,2), np.around(self.__time_remaining(dt,i,n),2)))
        self.__print("Finished sweep with total time {}".format(time.time()-t0))
        return outs
    
    def __time_remaining(self, dt: np.ndarray[int], i: int, n: int) -> int: #code
        """
        \ndt: case time array
        \ni: current case
        \nn: total cases
        \n
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

    def _invis_polar(self, a_array: np.ndarray, trim: str): #pre
        self.output_config(['t Alpha', 't CLtot', 't CDind', 't Elevator'])
        polar_dict = self.sweep(vc='a a', v_array=a_array, executes=[('vcv', 'd5 pm 0')])
        f = open('AVL/Output/Invis Polar.dat', 'w')
        f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format('Alpha', 'CLtot', 'CDind', 'Elevator'))
        f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format('-'*20, '-'*20, '-'*20, '-'*20))
        for k, val in polar_dict.items():
            f.write('{:>20}  {:>20}  {:>20}  {:>20}\n'.format(val['Alpha'], val['CLtot'], val['CDind'], val['Elevator']))

    def stall_prediction(self, Clx: float) -> float: #pre
        """
        \nStall prediction based on airfoil maximum Cl (Clx)
        \nChecks strip forces along wing. If strip Cl exceeds 95% Clx, the wing stalls
        """
        Clx = 0.95*Clx
        rdsub = RAMDaX(planename=self.__planename)
        rdsub.output_config(self.__output_list)
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
            rdsub._load()
            rdsub._oper()
            rdsub.vcv('a a {}'.format(a))
            rdsub._x()
            rdsub._save()
            rdsub._run()
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


class A_DaX():
    """
    \nParent class for design objects, since they all have a similar structure.
    \nContains functions common to all design objcets; does nothing on its own.
    """

    def __init__(self):
        pass

    "Modify an existing value in .config"
    def modify(self, **kwargs):
        keys = self.config.keys()
        for kw, val in kwargs.items():
            if kw in keys: 
                self.config[kw] = val
            else: 
                raise KeyError('Cannot modify {}: KeyError "{}"'.format(self, kw))
    
    "Only use this to append non-standard parameters to des objects"
    def assign(self, **kwargs):
        for kw, val in kwargs.items():
            self.config[kw] = val
        self.format += '1'
        self.kw += '1'

    "Used for __init__ assignments"
    def init_assign(self, kwargs):
        keys = self.config.keys()
        for kw, val in kwargs.items():
            if kw in keys: 
                self.config[kw] = val
            else: 
                raise KeyError('Cannot modify {}: KeyError "{}"'.format(self, kw))

    "Print label-value pairs to .avl files"
    def lv_pair(self, keys: str | list[str], cmt: bool = True):
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
    
    "Handles des object .config print order"
    def config_to_avl(self, f):
        keys = []
        for i, label in enumerate(self.config.keys()):
            keys.append(label)
            if self.format[i] == '1':
                if self.kw[i] == '1': cmt = False
                else: cmt = True
                f.write(self.lv_pair(keys=keys, cmt=cmt))
                keys = []
            else: continue

    "Print all children of a des object"
    def children(self, *args):
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

    "Print all parents of a des object"
    def parents(self, *args):
        print('Parent of {}:'.format(self))
        print('    {}'.format(self.parent))
    
class des_plane(A_DaX):
    """
    \nPlane/header component of .avl files
    \nParent: None
    \nChildren: des_surf
    """
    def __init__(self, **kwargs):
        "Direct AVL config"
        self.config = {'Name': 'Unnamed plane',
                       'Mach': 0.0,
                       'iYsym': (0), 'iZsym': int(0), 'Zsym': 0.0,
                       'Sref': 0.0, 'cref': 0.0, 'bref': 0.0,
                       'Xref': 0.0, 'Yref': 0.0, 'Zref': 0.0,
                       'CDp': 0.0}
        ".avl formatting"
        self.format = '110010010011'
        self.kw = '000000000000'
        "Track surfaces"
        self.surfs = []
        "Assignment"
        self.init_assign(kwargs)
        "Internal variables"
        self.reference = None
        "Rendering variables"
        self.surf_book = {}
        self.ctrl_book = {}

    "Custom printing"
    def __str__(self):
        return 'des_plane object "{}"'.format(self.config['Name'])
    
    "Begin writing the geometry to .avl file"
    def write_to_avl(self, filename: str):
        width = 95
        if self.reference:
            print('[RD] Using reference surface {} to update reference values'.format(self.reference))
            self.assign(Sref=self.reference.Sref, bref=self.reference.bref, cref=self.reference.Cref)
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
        print('done了')

    "WIP"
    def render(self):
        for sf in self.surfs:
            j = len(sf.secs) - 1
            pts = np.zeros((4*j+1,1))
            p = -1
            for i, sc in enumerate(sf.secs):
                a1 = int(-1/2 - p/2) 
                a2 = int(-1/2 + p/2)
                pts[2*i + a1] = sc.config['Xle']
                pts[2*i + a2] = sc.config['Xle'] + sc.config['Chord']
                pts[4*j-2*i + a1] = sc.config['Xle']
                pts[4*j-2*i + a2] = sc.config['Xle'] + sc.config['Chord']
                p *= -1
            self.surf_book[sf.config['Name']] = pts


class des_surf(A_DaX):
    """
    \nSurface component of .avl files
    \nParent: des_plane
    \nChildren: des_sec
    """
    count = 0
    def __init__(self, parent: des_plane, **kwargs):
        "Direct AVL config"
        self.config = {'Name': 'Unnamed surface',
                       'Nchordwise': 8, 'Cspace': 1.0, 'NSpanwise': 12, 'Sspace': 1.0,
                       'YDUPLICATE': None,
                       'TRANSLATE': (0.0,0.0,0.0),
                       'CDCL': None}
        ".avl file formatting"
        self.format = '10001111'
        self.kw = '00000111'
        "Coupling"
        self.parent = parent
        self.parent.surfs.append(self)
        self.secs = []
        "Assignment"
        self.init_assign(kwargs)
        "Internal variables"
        self.reference = []
        self.Sref = 0
        self.Cref = 0
        self.bref = 0
        self.count += 1

    "Custom printing"
    def __str__(self):
        return 'des_surf object "{}"'.format(self.config['Name'])

    "Calculate surface reference values - only used for custom reference"
    def assemble(self):
        Yle = [sc.config['Yle'] for sc in self.secs]
        Chord = np.array([sc.config['Chord'] for sc in self.secs])
        Yle = np.array(sorted(Yle))
        c = np.array([])
        y = np.array([])
        for i in range(len(Yle)-1):
            c = np.append(c, np.linspace(Chord[i], Chord[i+1], 101))
            y = np.append(y, np.linspace(Yle[i], Yle[i+1], 101))
        self.Sref = 2*np.trapz(c,y)
        self.Cref = 2*np.trapz(c**2, y)/self.Sref
        self.bref = 2*Yle[-1]
        # print(self.Sref, self.Cref, self.bref)

    "Sort surface sections into proper order"
    def order(self, orientation: int = 0):
        # 0 horizontal
        # 1 vertical
        Yle_collection = []
        for i, sc in enumerate(self.secs):
            Yle_collection.append([sc.config['Yle'], i, sc])
        Yle_sorted = sorted(Yle_collection, key=lambda x: x[0])
        for i, new_sc in enumerate(Yle_sorted):
            self.secs[i] = new_sc[2]

class des_sec(A_DaX):
    """
    \nSection component of .avl files
    \nParent: des_surf
    \nChildren: des_ctrl
    """
    def __init__(self, parent: des_surf, **kwargs):
        "Direct AVL config"
        self.config = {'Xle': 0.0, 'Yle': 0.0, 'Zle': 0.0, 'Chord': 1.0,
                       'Ainc': 0.0, 
                       'Nspanwise': 0, 'Sspace': 0,
                       'AFILE': None}
        ".avl file formatting"
        self.format = '00000011'
        self.kw = '00000001'
        "Coupling"
        self.parent = parent
        self.parent.secs.append(self)
        self.ctrls = []
        "Assignment"
        self.init_assign(kwargs)

    "Custom printing"
    def __str__(self):
        TRANSLATE = self.parent.config['TRANSLATE']
        if TRANSLATE is None:
            TRANSLATE = (0,0,0)
        else: pass
        global_pos = (self.config['Xle'] + TRANSLATE[0], 
                      self.config['Yle'] + TRANSLATE[1], 
                      self.config['Zle'] + TRANSLATE[2])
        return 'des_sec object located at (Xle, Yle, Zle) = ({:> 7g} {:> 7g} {:> 7g})'.format(*np.around(global_pos,3))
    
    "Create a new section using linear interpolation between two existing sections"
    def interp(s1, s2, span):
        val = []
        x1 = [s1.config['Yle'], s2.config['Yle']]
        for kw in ['Xle', 'Zle', 'Chord', 'Ainc']:
            x2 = [s1.config[kw], s2.config[kw]]
            val.append(np.interp(span, x1, x2))
        return des_sec(parent=s1.parent, Xle=val[0], Yle=span, Zle=val[1], Chord=val[2], Ainc=val[3], AFILE=s1.config['AFILE'])

class des_ctrl(A_DaX):
    """
    \nControl component of .avl files
    \nParent: des_sec
    \nChildren: None
    """
    def __init__(self, parent: des_sec,  **kwargs):
        "Direct AVL config"
        self.config = {'Cname': 'Unnamed control', 'Cgain': 1.0, 'Xhinge': 0.7, 'HingeVec': (0.0, 0.0, 0.0), 'SgnDup': 1.0}
        ".avl formatting"
        self.format = '00001'
        self.kw = '00000'
        "Coupling"
        self.parent = parent
        self.parent.ctrls.append(self)
        "Assignment"
        self.init_assign(kwargs)

    "Custom printing"
    def __str__(self):
        return 'des_ctrl object "{}"'.format(self.config['Cname'])

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
        \npt = np.array([[0.0,1.0,1.1,0.4,  0.0,1.0,1.1,0.4,  0.0],[0.0,0.0,2.4,2.4,0.0,0.0,-2.4,-2.4,0.0],[0.0,0.0,0.0,0.0, 0.0,0.0,0.0,0.0,0.0]])
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
        error_flags = ("Corrupted", "first!", "failed", "large", "zero-camber")
        error_book = {"Corrupted": "Error finding or loading .avl file. ",
                      "first!": "Flow execution not completed. Use _x() or check convergence.",
                      "failed": "Trim convergence failed. Check your variable-constraints.",
                      "large": "Resulting alpha is too large.",
                      "zero-camber": "Airfoil incorrectly loaded.",
                      "not recognized": "Unrecognized command passed to AVL."}
        log = open("Log.txt").read()
        for flag_term in error_flags:
            if re.search(flag_term, log):
                print("[Autodebug]", error_book[flag_term])
            else: pass
        