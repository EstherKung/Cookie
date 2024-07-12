from Code.Codebase import *

ad = ADaX(planename='dg500-extd.avl')
# ad.output_config(['t Alpha', 't e', 't CLtot', 't CDind', 's Xnp'])
# ad.inviscid_polar(np.arange(0,1.1,0.25), trim='d1 pm 0')
# ad.stall_prediction(Clx=1.4)
ad._load()
ad._oper()
# ad.use_run("Nu_Run_File.run", 1)
ad.vcv('d1 pm 0')
ad.vcv('a a 0')
ad._x()
ad._save()
ad.run(print_output=True)