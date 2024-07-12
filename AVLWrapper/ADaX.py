

"""""""Control Panel"""""""
Reload_Geometry = 1
Run_AVL         = 1
Autodebug       = 1
""""""""""""""""""""""""""





if Reload_Geometry: 
    import Code.__Aircraft_Geometry__ as __Aircraft_Geometry__
if Run_AVL:         
    import Code.__AVL_Commands__ as __AVL_Commands__
if Autodebug:       
    from Code.Codebase import tools
    tools._autodebug()
