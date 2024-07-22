#class Materials:
#    def __init__(this, name, vary = False):
#        this.name = name
# make into a class for future versions

def materials(des: str, relative = False):
    if not relative:
        match des:
            case "balsa": return 0.000022917375 #slugs/in^2
            case "foam": return 0.00001798666418 #slugs/in^3
            case "NLG": return 0.0084984484 #slugs
            case "MLG": return 0.01739492332 #slugs

            case "prop": return  0.002055654 #g // motor + esc + prop, we have 2
            case "battery": return  0.022269585 #slugs
            case "ESC": return 0.0045909606 #slugs
            case "motor": return 0.01305402 #slugs