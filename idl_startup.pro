if (!D.name EQ "X") then  device, GET_SCREEN_SIZE=screen_size $
		    else  screen_size = [ !D.x_size, !D.y_size ]

defsysv, "!DEVX", screen_size(0)-10
defsysv, "!DEVY", screen_size(1)-24
defsysv, "!DEBUG", 0

!PATH = "~/idl/vlib/widgets:"	+ $
	"~/idl/vlib/math:"	+ $
	"~/idl/vlib/image:"	+ $
	"~/idl/vlib/color:"	+ $
	"~/idl/vlib/window:"	+ $
	"~/idl/vlib/graphics:"	+ $
	"~/idl/vlib/struct:"	+ $
	"~/idl/vlib/misc:"	+ $
	"~/idl/vlib/io:"	+ $
	"~/idl/vlib/data:"	+ $
	"~/idl/vlibm/mgep:"	+ $
	"~/idl/vlibm/rtfit:"	+ $
	"~/idl/vlibm/models:"	+ $
	"~/idl/vlibm/deconv:"	+ !PATH

if (!D.name EQ "X") then !PATH = "~/idl/vlib/xdoc:" + !PATH

const = {CONSTANTS_MKS, $
         alpha: 7.2973525698e-3,$;fine structure constant
         as2sr: 2.35044305391e-11,$; arcsec^2 --> sr conversion factor
         AU:149597870700.,$      ;AU, in m (exact)
         c: 299792458. ,$        ;vacuum speed of light, m/s
         dtor: 0.0174532925199,$ ;pi/180 deg-->rad conversion factor
         d2tosr: 3.046174198e-4,$;(pi/180)^2, deg^2-->sr conv. factor
         e: 2.7182818284590452 ,$;Euler's number                       
         ev: 1.602176565e-19,$   ;elementary charge e, 1 electron volt, C
         eps0: 8.854187817e-12 ,$;electric vacuum permittivity, F/m
         F: 96485.3365,$         ;Faraday constant NAe, C/mol
         G: 6.67384e-11,$        ;Gravitation constant, m3/kg/s2
         gn: 9.80665,$           ;Earth standard gravity, m/s2
         h: 6.62606957e-34,$     ;Planck constant, J s
         hbar: 1.054571726e-34,$ ;h/(2pi), J s
         k: 1.3806488e-23,$      ;Boltzmann constant R/NA, J/K
         Lsun: 3.828e+26 ,$      ;solar luminosity, Watts
         Msun: 1.98855e+30 ,$    ;solar mass, kg.
         Rsun: 6.5966e+8 ,$      ;solar radius, m.
         me: 9.10938291e-31,$    ;electron mass, kg
         mn: 1.674927351e-27,$   ;neutron mass, kg
         mp: 1.672621777e-27,$   ;proton mass, kg
	 mH: 1.673534e-27,$	 ;mass of neutral H atom, kg
         mu0: 12.566370614e-7,$  ;magnetic vacuum permeability, N/A2
         Na: 6.02214129e23,$     ;Avogadro constant NA, mol-1
         pc: 3.08567758065e+16 ,$;parsec, in m
         pctoAU:206265.,$        ;pc to AU conversion factor
         phi: 1.6180339887498948,$ ;golden ratio
         pi: 3.1415926535897932,$  ;Pi 
         R: 8.3144621,$          ;molar gas constant, J/mol/K
         R_earth: 6370997.0,$    ;Earth radius (spherical), m
         re: 2.8179403267e-15,$  ;classical electron radius, m
         ryd: 10973731.568539,$  ;Rydberg constant Rinf, m-1
         sigma: 5.670373e-8,$    ;Stefan-Boltzmann constant, W/m2/K4
         u: 1.660538921e-27}     ;unified atomic mass unit, kg 

cv = {	CONSTANTS_CGS,	$
	h: 6.6260755e-27 ,$ ;Planck's constant, erg-sec.
	sb: 5.67051e-5	 ,$ ;Stephan-Boltzman constant, erg/cm^2/deg^4.
	k: 1.380658e-16	 ,$ ;Boltzman constant, erg/Kelvin.
	wtmax: 2897.3	 ,$ ;micron-deg, wavelen * max( T ) of B-B spectrum.
	c: 2.99792458e10 ,$ ;cm/sec, speed of light in vacuum.
	ev: 1.602177e-12 ,$ ;ergs/eV.
	g: 6.67e-8	 ,$ ;Gravitation constant, cm^3/gm/sec^2.
	mH: 1.6733e-24	 ,$ ;mass of Hydrogen, gm.
	mp: 1.672623e-24 ,$ ;mass of proton, gm.
	me: 9.10939e-28	 ,$ ;mass of electron, gm.
	mn: 1.674929e-24 ,$ ;mass of neutron, gm.
        amu: 1.66054e-24 ,$ ;atomic mass unit (m_carbon/12), gm.
	Lsun: 3.826e33	 ,$ ;solar luminosity, ergs/sec
	Msun: 1.989e33	 ,$ ;solar mass, gm.
	Rsun: 6.596e10	 ,$ ;solar radius, cm.
	a0: 5.292e-9	 ,$ ;Bohr radius, cm.
        fsc: 7.297353e-3 ,$;fine structure constant. 
	txsec: 6.65e-25  ,$ ;Thomson cross section, cm^2.
        Ryd: 2.17987e-11 ,$ ;Rydberg constant, ergs.
	pc: 3.086e18	 ,$ ;parsec, cm.
	AU: 1.496e13	 }  ;sun to earth distance, cm.

csolgal = {CONSTANTS_SOLGAL, $
        G: 4.302e-3    ,$ ;Gravitation constant, (pc*(km/s)^2)/Msun
        AU: 4.8481e-6  ,$ ;sun-earth distance, pc
        LY: 0.306601   ,$ ;LY, in pc
        vsol: 255      ,$ ;circular orbit speed of the sun, km/s
        vlsr: 220      ,$ ;local standard of rest velocity, km/s*
        Mmw: 8.5e+11   ,$ ;Milky May mass, Msun**
        R0: 7.98e+3    ,$ ;Solar Circle radius, pc
        Rgal: 15.0e+3  ,$ ;Milky Way radius, pc (assumes no 23-27.5 kpc ring)
        H0: 67.8       ,$ ;Hubble constant, km/s/Mpc.
        tH: 1.44e+10   ,$ ;Hubble time (1/H0), yrs.
        LH: 4.550e+9   } ;Hubble length (cH0), pc.
; *(accepted values range from 202 to 241 km/s)
;**(accepted values range from 5.8 to 15 e+11 Msun)

defsysv, "!CV", cv, 1	;define as read-only sys-var.
defsysv, "!CONST", const, 1	;define as read-only sys-var.
defsysv, "!CSOLGAL", csolgal, 1	;define as read-only sys-var.

@~/idl/setup_user	;invoke a user defined global setup.

@setup	;invoke whatever setup is present in current directory,
	; otherwise, ignore the error message about file not found.
