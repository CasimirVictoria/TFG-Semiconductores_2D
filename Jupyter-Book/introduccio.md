---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.10.3
kernelspec:
  display_name: SageMath 9.2
  language: sage
  name: sagemath
---

# Introducció

```{note}
Per ara sols són algunes idees a desenvolupar, fins que no acabe els càlculs no vaig a posar-me a redactar la memòria, però mentres tant vaig provant les funcionalitats combinades de sagemath, jupyter, jupytext i jupyterbook.
```

## Model de Born i Von Karman
- Explicar el **model de  Born y  Von Karman (1921)**, tal y com ve en {% cite brueesch82_phonon %} (pàgina 4).


## Aproximació adiabàtica
Desprès explicar que fent ús de l'**aproximació adiabàtica**, podem, sempre que siga vàlida tal aproximació -explicar quan ho és- considerar els nuclis iònics y els electrons de valència com constituents del sòlid independents, de manera que podem escriure l'energia potencial (o el Hamiltonià) como una suma de les distintes contribucions.
  
## Aproximació harmònica
Passar a descriure l'**aproximació harmònica**: escriure el potencial, les equacions de moviment de Lagrange, y passar d'estes equacions de moviment en l'espai de posiciones al problema de valors propis de la matriu dinàmica en l'espai de vectores de onda (o moments).

Les vibracions reticulars estan regides per les forces que experimenten els àtoms quan es desplacen de la seua posició d'equilibri. La primera hipòtesi és que cada àtom té una posició d'equilibri en el cristall, y considerarem que estos àtoms vibren amb una amplitud menuda (en comparació amb la distància interatòmica) al voltant d'aquesta, de manera que el sòlid es troba en estats que corresponen al que macroscòpicament es coneix com *la regió de comportament el·lastic lineal*, on es verifica la llei de Hooke.

Podem per tant aproximar l'energia potencial de interacció pel terme armònic del seu desenvolupament en sèrie de potencies del desplaçament ....

- Explicar el concepte de fonó i la seua importància: Un fonó és un mode de vibració quantitzat que té lloc en una xarxa cristallina. Els fonons tenen una gran importància en moltes propietats físiques dels sòlids. Los fonons son una version mecano-quàntica dels modes normals de vibració de la mecànica clàssica, on cada par de la xarxa osci\l.la amb la mateixa freqüència. Aquests modes normals són importants perquè qualsevol moviment vibracional de la xarxa pot descriure's com una superposició de modes normals de distinta freqüència, en este sentit són la base de les vibracions de la xarxa.

- Explicar també que una vegada descrit el cristall en l'espai de posicions tenim que descriure la xarxa cristallina en l'espai de moments ...

- Una vegada conclosa la introducció teòrica general, descriure el cristall de BN monocapa (xarxa cristallina, base, xarxa en l'espai $2\pi\text{-recíproc}$.

- Passar a construir la matriu dinàmica a partir del tensor de constants de for\c{c}a (cite Balkanski 2000) en coordenades cartesianes, tenint en compte que degut a la simetria del cristal la matriz de constantes de fuerza (y su transformada de Fourier, la matriz dinámica), té que complir certes relacions (forma de la matriz, elementos nulos, relación entre elementos, ...), cite aizawa90 bond soften monol graph formed (en el apèndix). 


## Xarxa recíproca. Primera zona de Brillouin
Donada una terna de vectors base primitius, $p_i$, de la xarxa cristallina en l'espai de posicion, la condició:

```{math}
:label: rec-1
\vec p_i\cdot\vec p_j^{*}=\delta_{ij}
```

on $\delta_{ij}$ és la delta de Kronecker, pot considerar-se un sistema d'equacions que defineix una altra terna de vectors $p_j^*$. La solució de l'equació {eq}`rec-1` és

```{math}
:label: rec-2
\vec p_1^*=\frac{\vec p_2\times\vec p_3}{\vec p_1\cdot(\vec p_2\times\vec p_3)},\quad \vec p_2^*=\frac{\vec p_3\times\vec p_1}{\vec p_1\cdot(\vec p_2\times\vec p_3)}, \quad \vec p_3^*=\frac{\vec p_1\times\vec p_2}{\vec p_1\cdot(\vec p_2\times\vec p_3)}
```

Els *vectors recíprocs* $\vec p_j^{*}$ són els vectors d'una altra xarxa cristallina coneguda com *xarxa recíproca*. Els  vectors recíprocs, i els vectors de translació cristallina recíprocs, tenen dimensions de inversa de longitud i es representen en l'*espaci recíproc* o de números d'ona. **Les xarxes cristallines real y recíproca són dues descripcions equivalents del matexi sistema físic: el sòlido cristallí** que s'està estudiant.


Tenim que tindre present que l'espai $2\pi$-recíproc és el fonamental en l'estudi de la Física de l'Estado Sòlid ja que els estats de les partícules i les interaccions físiques de interès es descriuen en l'espai de vectors d'ona.

Podem interpretar que 

```{math}
:label: rec-3
\vec K_{h_1h_2h_3}=2\pi\left(h_1\vec a_1^{*}+h_2\vec a_2^{*}+h_3\vec a_3^{*}\right)
```

són els vectors de translació cristallina que defineixen una xarxa cristallina en aquest espai $2\pi$-recíproc (sols es diferencia de l'espai recíproc en un factor d'escala $2\pi$). En termes físics, l'espai $2\pi$-recíproc és l'/espai de vectors d'ona/ $\vec k$, i a falta d'un factors d'escala $\hbar$, coincideix amb l'espai de moments  $\vec p=\hbar\vec k$. 

### Cel·les de Wigner-Seitz (WS)
Para descrire la xarxa 2$\pi$-recíproca, emprarem el criteri de *Wigner-Seitz*. Les *cel·les de Wigner-Seitz* estan centrades en un nuc de la xarxa i es defineixen com la regió més pròxima a un nuc (el del centre de la cel·la) que a qualsevol altre. Per determinar la seua forma, partim d'un nuc qualsevol com a origen, construim els segments que uneixen este nuc amb els seus veïns i es tracen els plans que bisecten cadascun d'estos segments: la cel·la de WS és la cel·la de menor volum al voltant de l'origen que està delimitada per estos plans (rectes en el caso d'una xarxa bidimensional).

Notem que en l'espai de $\vec k$ emprarem cel·les unitat de WS mentre que en l'espai de posicions sempre emprem cel·les unitat de Bravais.
La cel·la de WS de la xarxa $2\pi$-recíproca es coneix com **primera zona de Brillouin** (ZB)

## Matriu dinàmica

La matriu dinàmica és la magnitud central de la dinàmica reticular: les freqüències dels fonons es calculen a partir dels valors propis de la matriu dinàmica:

```{math}
\sum_{\alpha\prime}D_{\alpha\alpha\prime}\cdot\vec e_{\alpha\prime}(\vec q)=\omega^{2}\vec e_{\alpha}(\vec q)
```

Per tant, les freqüències $\omega$ com funció del vector d'ones $\vec q$ del fonó són solució de l'equació secular:

```{math}
\det\left|\frac{1}{\sqrt{M_\alpha M_{\alpha\prime}}}D^{ij}_{\alpha\alpha\prime}\left(\vec q\right)-\omega^2\left(\vec q\right)\right|=0 
```

on $M_{\alpha}$ es la massa de l'àtom $\alpha$ y la matriu dinàmica ve definida por:

```{math}
:label: Matriu_Dinàmica
D_{\alpha,\alpha\prime}^{i,j}=\frac{\partial^2 E}{\partial u_{\alpha}^{*i}(\vec q)\partial u_{\alpha\prime}^{j}(\vec q)}
```

on $u_{\alpha}^{i}$ representa el desplaçament de l'àtom $\alpha$ en la direcció $i$.

La segona derivada de l'energia de l'equació {eq}`Matriu_Dinàmica` correspon al canvi en la força que actua sobre l'àtom $\alpha\prime$ en la direcció $j$ quan es desplaça l'àtom $\alpha$ en la direcció $i$

```{math}
D_{\alpha\alpha\prime}^{ij}(\vec q)=\frac{\partial}{\partial u^{*\alpha}_{i}}F^{j}_{\alpha\prime}(\vec q)
```
