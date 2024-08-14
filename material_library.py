#class Materials:
#    def __init__(this, name, vary = False):
#        this.name = name
# make into a class for future versions

def materials(des: str, relative = False):
    if not relative:
        match des:
            case "balsa": return 0.000022917375 #slugs/in^2
            case "foam": return 0.00001798666418 #slugs/in^3
            case "NLG": return 4.91e-3 #slugs
            case "MLG": return 0.006038601994 #slugs

            case "prop": return  0.0031469472 #slugs
            case "battery": return  0.022269585 #slugs
            case "ESC": return 0.0045909606 #slugs
            case "motor": return 0.01305402 #slugs
            case "Rx batt": return 0.009593052 #slugs
            case "Rx": return 0.0009593052 #slugs

            case "payload": return 0.0 #slugs, pretend value