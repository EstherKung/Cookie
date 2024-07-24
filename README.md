## Limitations of this program

Masses that are not (yet) accounted for

- servos, wire
- adhesives
- connection between empennage and fuselage

Improvements to be made on estimations: 

- Wing— currently uses just foam
- Fuselage— Dean has it work in progress

## Configuration

### Wing

No quarter chord sweep for the main wing

only one section, constant taper

taper ratio fixed at $\lambda=0.4$

### Landing Gear

Assumes NLG is placed at the start of the main fuselage

MLG placement 10% of fuselage length behind the CG— does not include empennage length.

used the masses provided by David in [Landing Gear](https://www.notion.so/Landing-Gear-d85f9a39079147799af1fb694bdc6855?pvs=21) 

should have another option to implement the percent estimates at more initial stages of design

- landing gear is 10% total weight
    - 30% at the NLG
    - 70% at the MLG

### Fuselage

Specifically made for a fuselage configuration with three sections: 

- tapered nose section
- constant fuselage section
- tapered aft section

Assumes the main fuselage section ends where the wing root chord ends

(inherently) Assume top surfaces of the fuselage is in plane

### Empennage

Assumes the leading edge of the empennage is attached to the end of the fuselage section

Assume hstab is in plane with the wing

## Code Structure

The following files here contains class definitions to initiate the plane from its components.

`plane_def` 

has a general parent class, Components

```python
"""
Attributes: mass, cg, location, and geometry parameters

Re: CG/Mass Buildup
materials are called out from material_library.py
depending on type, the mass is calculated from material density

Re: Drag Tally
default components do not have S_wet and FF is set to 1
"""

"has methods"
m() #returns mass of component
"Re: drag tally"
wetted_area() #defines the wetted area attribute of the specified part. 
FF_calc() #calculates the form factor of the wing, or the fuselage, depending on designation.
drag_visc() #returns the c_d of the component
"for plotting"
geometry() #returns a rectangle of the part (currently only used by the empennage)
```

there are also components that have more specific geometry in child classes 

`Wing`

made such that all parameters of planform can be defined, along with the four angles

utilises methods in `wing_def`

mass is calculated from the cross section profile across the span, integrated with trapezoidal rule

cg is derived from geometry centroid of the wing

thickness is *thickness to chord*

overwrites `geometry` to plot taper

`Fuselage`

made to interface with `fuse_def` file made by Dean. 

specifies the fuselage section as `nose`, `main`, and `aft`

define the length of this section, and their cross section in a pair of list

overwrites method `geometry()` to plot the tapered section

`Landing Gear`

overwrites method geometry() to plot landing gear as a little square. 

`wing_def`

has methods used to aid the component definition of the wing. 

these here take argument of the chosen airfoil dat file

`csec(afile)`: finds the cross section area (for unit chord) of the airfoil

`centroid(afile)`: finds the geometric centroid (for unit chord) of the airfoil

`thickness(afile)`: finds the thickness (for unit chord) of the airfoil

`camber(afile)`: [INCOMPLETE] intend to calculate the camber of the airfoil

and this is from Cameron homework 

`planform()`: spits out all the planform parameters. takes arguments in a list of order `[AR, b, S, taper_ratio, c_r, c_t]`

`xCG_calc` and `nearest()` are sub functions to aid calculation, not used directly to define the wing.

`fuse_def`

that’s Dean’s

`material_library`

all the materials and their respective density. 

`sauce`

where we run the sauce. The routine goes 

1. pick a point on the constraint diagram
2. Mass Convergence Study
3. CG placement Study (vary wing location)
4. Stability Analysis in AVL (vary empennage length)