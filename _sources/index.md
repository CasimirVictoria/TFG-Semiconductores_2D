# Sobre este treball

## Notes sobre la versió en línia


```{note}
La versió en línia continuaré actualitzant-la, per provar i incloure funcionalitats de *jupyterbook* (comentaris emprant *hypothesis*,  etc.)
```

## Resum 

Els materials bidimensionals (2D) com el grafè són de gran interès tant per les seues propietats físiques exclusives com per les seues aplicacions potencials. L'estudi de la dinàmica de la xarxa cristal·lina (*fonons*) d'estos materials és un requisit previ per entendre la seua estabilitat estructural i propietats tèrmiques, així com les seues propietats de transport i òptiques. 
 
 En este treball es calculen els modes vibracionals  d'un  material semiconductor 2D model, el nitrur de bor monocapa. Emprant un model clàssic i senzill, el model de constants de força, s'obtenen les expressions analítiques per als modes de vibració als punts crítics de la zona de Brillouin. Evidentment, un model tan senzill com l'emprat, on sols es consideren les interaccions amb uns pocs veïns i on s'assumeix a més una forma diagonal per al tensor de constants de força, té moltes limitacions, però ajustant les constants de força emprant la relació de dispersió calculada a partir de primers principis mitjançant el mètode de *Density Functional Perturbation Theory* observem que, almenys, obtenim resultats acceptables per a baixes freqüències.

## Abstract
Two-dimensional materials (2D) like graphene are of great interest for both their unique physical properties and their potential applications. The study of the dynamics of the crystal lattice (*phonons*) of these materials is a prerequisite for understanding their structural stability and thermal properties, as well as their transport and optical properties.
 
 In this work the vibrational modes of a 2D model semiconductor material, the monolayer boron nitride, are calculated. Using a classical and simple model, the model of force constants, analytical expressions for the modes at the critical points of the Brillouin zone are obtained. Obviously, a model as simple, where only interactions with a few neighbors are considered and where a diagonal shape is assumed for the force constant tensor, has many limitations, however by fitting the force-constants to the ab initio  dispersion relation, calculated using the method of *Density Functional Perturbation Theory* we observe that, at least, we obtain acceptable results for low frequencies.
 

  
## Eines emprades

Els càlculs s'han realitzat en el llenguatge de programació [python](https://www.python.org/){cite:p}`4160250`, emprant el software matemàtic [sagemath](https://www.sagemath.org/) {cite:p}`sagemath`. Sagemath està construït sobre python i facilita un accés unificat a biblioteques orientades a càlcul numeric com [numpy](https://numpy.org/) {cite:p}`harris2020array` , [scipy](https://scipy.org/) {cite:p}`2020SciPy-NMeth`, biblioteques orientades a càlcul simbòlic com [sympy](https://www.sympy.org/) i el sistema d'àlgebra computacional [maxima](https://maxima.sourceforge.io/), entre molts altres. Si es necessita emprar python o qualsevol de les seues biblioteques es poden emprar directament també; en este treball, per exemple, he emprat també la biblioteca [pandas](https://pandas.pydata.org).

Com a entorn de desenvolupament s'ha emprat [jupyterlab](https://jupyter.org/) i [git](https://git-scm.com/) com sistema de control de versions distribuït. 

El projecte, incloent els *jupyter notebooks* amb els càlculs, es troba accessible al repositori públic de github [TFG-Semiconductores\_2D](https://github.com/CasimirVictoria/TFG-Semiconductores_2D). Així mateix, una versió en línia, orientada a mostrar els càlculs realitzats, es troba publicada en [github pages](https://casimirvictoria.github.io/TFG-Semiconductores_2D/index.html); la creació i publicació d'esta versió online s'ha realitzat automàticament a partir dels *jupyter notebooks* emprant [jupyterbook](https://jupyterbook.org/intro.html).

Jupyter/Jupyterlab ha demostrat ser un entorn fantàstic (i molt popular) per a programar en python (i altres llenguatges), però al no treballar sobre fitxers de text pla (en este cas fitxers amb codi python), dificulta emprar eines de control de versions de codi com git. Per aquesta raó s’ha emprat també el plugin [jupytext](https://jupytext.readthedocs.io/en/latest/) per a Jupyter, que permet tindre sincronitzat el notebook amb un fitxer en text pla (amb el codi en python, markdown, etc.). En este projecte he decidit emprar el format [myst markdown](https://jupyterbook.org/content/myst.html), de manera que des de jupyterlab obric un fitxer markdown i la extensió jupytext sincronitza un fitxer en format *.ipynb* amb este fitxer, què és amb el qual realment treballa jupyterlab. Este *notebook* l'exporte després a un fitxer python que carregue posteriorment al fitxer .tex on he escrit la memòria, emprant el paquet [sagetex](https://ctan.org/pkg/sagetex), de manera que puc incloure tots els càlculs realitzats (incloent expressions analítiques, no sols numèriques i gràfics) de manera automàtica en la memòria sense tindre que escriure-les manualment, aprofitant l'excel·lent suport que té sagemath de $\LaTeX$ . 

Al repositori de github, dins de la carpeta **Memoria** es troba un petit *script* en *bash* que automatitza la generació d'esta memòria si realitze algun canvi sobre el notebook amb els càlculs.

Notem que tot els software emprat per a la realització del treball és de codi obert i gratuït.
