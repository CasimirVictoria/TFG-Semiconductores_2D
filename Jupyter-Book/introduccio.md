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

# Introducció teòrica

En este treball s'empra el model de constants de força (*FCM* pel seu acrònim en anglés, *force-constant method*), per descriure la dinàmica de la xarxa cristal·lina d'un material bidimensional, el **nitrur de bor** monocapa, *BN*. La metodología emprada per realitzar els càlculs es basa sobretot en els articles recomanats pels tutors, en particular {cite:p}`wirtz04_phonon_disper_graph_revis` i {cite:p}`falkovsky08_symmet_const_phonon_disper_graph`, tot i que en estos articles es tracta el grafé, i en el present treball estudiem el nitrur de bor.

En este model la dinàmica dels àtoms es descriu considerant que cada àtom interacciona amb els del seu entorn fins a un determinat nombre de veïns (els veïns es classifiquen segons la seua distància), caracteritzant les interaccions amb els veïns mitjançant un tensor de constants de forces, una manera visual de representar el model seria considerar que els àtoms estan conectats per molls. Altres mètodes, com els càlculs realitzats a partir de primers principis (*ab initio*) emprant la teoria del funcional densitat (*DFT*, pel seu acrònim en anglés *density-functional density*) o el mètode de camp de força de valència (*VFF*, *valence force field*) necessiten temps de càlcul molt més llargs (i donen resultats molt més precissos, clar). El mètode de constants de força empra un conjunt reduït de paràmetres que poden ajustar-se a les mesures experimentals o bé, com s'ha fet en este treball, a la relació de dispersió calculada per primers principis, pel mètode de *Density Functional Perturbation Theory* em este treball.  Notem que encara que estem tractant amb un mètode senzill, este ha demostrat que pot proporcionar resultats fiables {cite:p}`wirtz04_phonon_disper_graph_revis` .

Abans d'explicar amb detall el mètode emprat en este treball repassem alguns conceptes i aclarim la notació emprada.

## Xarxa cristal·lina

Sabem que els àtoms que constitueixen un sòlid cristal·lí estan distribuïts d'una manera regular en l'espai, y per descriure esta distribució regular s'empra el concepte matemàtic de *xarxa cristal·lina*, introduït per *A. Bravais* en 1845: una xarxa cristal·lina es defineix com una distribució discreta i regular de punts (*nucs*) que té sempre la mateixa aparença independentment de quin punt escollim com origen de coordenades.

En el cas que anem a tractar, un cristall bidimensional, la xarxa cristal·lina bidimensional pot generar-se a partir de dos vectors base, de manera que els vectors de posició dels *nucs de la xarxa*, els vectors de translació cristal·lina,  tenen la forma:

\begin{equation}
 \vec R_\vec l=\vec R_{l_1 l_2}=l_1 \vec a_1 +l_2\vec a_2
\end{equation}

on $\vec l$ es denomina índex del nuc (notem que en la literatura sobre el tema sol emprar-se com a índex la lletra $n$, però crec que un índex vectorial especifica millor els nucs). Si tots els nucs de la xarxa tenen índexs enters llavors els vectors base s'anomenen *primitius*.

Els cristalls o *estructures cristal·lines* són distribucions regulars d'àtoms en l'espai de posicions i els podem descriure associant a cadascun dels nucs d'una xarxa cristal·lina en l'espai de posicions un conjunt d'àtoms o base (matemàticament és el producte de *convolució* de la base y la xarxa de nucs). En el cas del *BN* tractem amb una base diatòmica, formada per un nucli de bor, $B$, i un de nitrogen, $N$. 

En dues dimensions, els paral·lelograms amb els quals podem omplir completament l'espai per translació cristal·lina i que contenen almenys un nuc de la xarxa es coneixen com cel·les unitat. La cel·la unitat més sencilla és la que té per costats el vectors base, i es coneixen com *cel·les de Bravais*.

### Xarxa recíproca. Primera zona de Brillouin
Donat un conjunt de vectors base, $p_i$, de la xarxa cristal·lina en l'espai de posicions, la condició:

\begin{equation}
\label{eq:rec1}
\vec p_i\cdot\vec p_j^{*}=\delta_{ij}
\end{equation}

on $\delta_{ij}$ és la delta de Kronecker, defineix un altre conjunt de vectors $p_j^*$, coneguts com vectors recíprocs, i què són els vectors base que defineixen una altra xarxa coneguda com *xarxa recíproca*. Els  vectors base recíprocs, i els vectors de translació cristal·lina recíprocs, tenen dimensions de inversa de longitud i es representen en l'*espai recíproc* o de nombres d'ona. **Les xarxes cristal·lines real y recíproca són dues descripcions equivalents** del mateix sistema físic: el sòlid cristal·lí que s'està estudiant.

Podem interpretar que
\begin{equation}
\label{eq:rec3}
\vec K_{h_1h_2}=2\pi\left(h_1\vec a_1^{*}+h_2\vec a_2^{*}\right)
\end{equation}

són els vectors de translació cristal·lina que defineixen una xarxa cristal·lina en l'espai $2\pi$-recíproc (sols es diferencia de l'espai recíproc en un factor d'escala $2\pi$). En termes físics, l'espai $2\pi$-recíproc és l'*espai de vectors d'ona* $\vec k$, i a falta d'un factors d'escala $\hbar$, coincideix amb l'espai de moments  $\vec p=\hbar\vec k$.

Cal tenir present que l'espai $2\pi$-recíproc és el fonamental en l'estudi dels sòlids cristal·lins, ja que els estats de les partícules i les interaccions físiques de interès es descriuen en l'espai de vectors d'ona, $\vec k$.


<p><b>Cel·les de Wigner-Seitz (WS) y primera zona de Brillouin (ZB)</b></p>

Per descriure la xarxa $2\pi$-recíproca, s'empra el criteri de *Wigner-Seitz*. Les *cel·les de Wigner-Seitz* estan centrades en un nuc de la xarxa i es defineixen com la regió més pròxima a un nuc (el del centre de la cel·la) que a qualsevol altre. Per determinar la seua forma, partim d'un nuc qualsevol com a origen, construïm els segments que uneixen este nuc amb els seus veïns i es tracen els plans que bisecten cadascun d'estos segments: la cel·la de *WS* és la cel·la de menor volum al voltant de l'origen que està delimitada per estos plans (rectes en el cas d'una xarxa bidimensional).

Notem que en l'espai de $\vec k$ s'empren cel·les unitat de *WS* mentre que en l'espai de posicions sempre emprem cel·les unitat de Bravais.
La cel·la de WS de la xarxa $2\pi$-recíproca es coneix com **primera zona de Brillouin**


## Vibracions atòmiques en cristalls
Passem ara a descriure el model emprat per descriure les vibracions dels àtoms del cristall.

### Model de BORN i VON KARMAN}

Els àtoms constituents d'un sòlid cristal·lí no estan immòbils sinó que vibren al voltant de la seua posició de equilibri.

En 1912, Born y Von Karman {cite:p}`Born:1912:SRG` introduïren un model per explicar la dinàmica cristal·lina, on la idea fonamental és que cada mode normal té l'energia de un oscil·lador de Planck.

Assumiren que els àtoms es troben disposats en una matriu tridimensional periòdica, i que la força sobre un àtom depèn no del desplaçament d'este respecte de la seua posició d'equilibri, sinó del desplaçament relatiu d'este àtom respecte als seus veïns. Introduïren també les condicions de contorn periòdiques que simplifiquen considerablement els càlculs. 

La dinàmica del sistema es descriu de manera senzilla, no en termes de les vibracions de àtoms individuals, sinó en termes de moviments col·lectius en forma de ones viatgeres anomenades vibracions cristal·lines (*lattice vibrations*) per Born. Cada vibració cristal·lina es caracteritza per una freqüència y un vector de ona.
La vibració cristal·lina quantitzada, o mode normal, s'anomena fonó per analogia amb el fotó.

Un fonó és un mode normal de vibració quantitzat que té lloc en una xarxa cristal·lina.  Aquests modes normals són importants perquè qualsevol moviment vibracional de la xarxa pot descriure's com una superposició de modes normals de distinta freqüència, en este sentit són la base de les vibracions de la xarxa. Els fonons tenen una gran importància en moltes propietats físiques dels sòlids. 

## Aproximació adiabàtica

Moltes propietats físiques dels sòlids poden classificar-se bé com electròniques o bé com vibracionals, segons estiguen determinades pels electrons (de valència) o per les vibracions dels àtoms (nuclis iònics): normalment considerem els nuclis i els electrons com a constituent independents del sòlid, ja que les masses dels electrons i dels nuclis són tan diferents que el moviment dels nuclis és molt més lent que el dels electrons. L *aproximació adiabàtica*, introduïda per Born i Oppenheimer {cite:p}`ANDP:ANDP19273892002` estableix que els electrons responen de manera pràcticament instantània als desplaçaments atòmics, de manera que el seu estat ve determinat per los posicions atómiques instantànies, mentre que els àtoms no poden respondre a les ràpides variacions espacials dels electrons: diguem que els electrons segueixen el moviment iònic adiabàticament.

Notem, però, que estem tractant amb una primera aproximació al problema y algunes propietats físiques venen determinades per la interacció entre els electrons y les vibracions atòmiques.


## Aproximació harmònica

Les vibracions cristal·lines estan regides per les forces que experimenten els àtoms quan es desplacen de la seua posició d'equilibri. La primera hipòtesi és que cada àtom té una posició d'equilibri en el cristall, que denotem per $\vec r^{(0)}_{\vec l,\vec\alpha}=\vec R_\vec l+\vec d_{\alpha}\quad $, y considerarem que estos àtoms vibren al voltant d'esta posició d'equilibri, $\vec r_{\vec l,\vec\alpha}=\vec r^{(0)}_{\vec l,\vec\alpha}+\vec u_{\vec l,\alpha}(t)\quad$, amb una amplitud menuda (en comparació amb la distància interatòmica) al voltant d'aquesta,   de manera que el sòlid es troba en estats que corresponen al que macroscòpicament es coneix com \textit{la regió de comportament elàstic lineal}, on es verifica la llei de Hooke.

Podem, per tant, aproximar l'energia potencial de interacció pel terme harmònic del seu desenvolupament en sèrie de potencies del despla\c{c}ament:

\begin{equation}
 V=\sum_{\vec l',\alpha'}\sum_{\vec l,\alpha}\frac{1}{2}\vec u_{\vec l'\alpha'}\cdot\underbrace{\mathbf\Phi^{\alpha',\alpha}(\vec R_{\vec l}-\vec R_{\vec l'})}_\text{matriu de constants de forces}\cdot\vec u_{\vec l,\alpha}
\end{equation}

Les equacions de moviment en l'aproximació harmònica s'escriuen en la coneguda forma:

```{math}
:label: movharmonic
 M_{\alpha'}\vec{\ddot{u}}_{\vec l,\alpha'}(t)=-\sum_{\vec l',\alpha}\mathbf\Phi^{\alpha',\alpha}\left(\vec R_\vec l-\vec R_{\vec l'}\right)\cdot\vec u_{\vec l,\alpha}(t)
```

 Aquesta equació representa un sistema  d'oscil·ladors harmònics acoblats, on $\alpha$ i $\alpha'$ fan referència al àtoms de la base considerats, l'índex vectorial $\vec l$ (que sol aparèixer en la literatura com $n$) índica el nuc considerat i $\vec R_\vec l$ és el vector de translació cristal·lina.

### Matriu dinàmica
Realitzant la transformada de Fourier de l'equació {eq}`movharmonic`, passem al conegut problema d'autovalors:

\begin{equation}
\sum_{\alpha}\mathbf D^{\alpha'\alpha}(\vec q)\cdot\vec e_{\alpha}(\vec q)={\omega'}^{2}\vec e_{\alpha'}(\vec q)
\end{equation}   

La matriu dinàmica és la magnitud central de la dinàmica reticular: les freqüències dels fonons es calculen a partir dels valors propis de la matriu dinàmica.

Per tant, les freqüències $\omega$ com funció del vector d'ones $\vec q$ del fonó són solució de l'equació secular:

\begin{equation}
\label{eq:secular}
\det\left|\frac{1}{\sqrt{M_{\alpha'} M_{\alpha}}}\vec D^{\alpha'\alpha}\left(\vec q\right)-{\omega'}^2\left(\vec q\right)\right|=0 
\end{equation}

on $M_{\alpha}$ es la massa de l'àtom $\alpha$ y la matriu dinàmica ve definida per:

```{math}
:label: Matriu_Dinàmica
\mathbf D^{\alpha',\alpha}_{i,j}=\frac{\partial^2 E}{\partial u^{\alpha'*}_i(\vec q)\partial u^{\alpha}_{j}(\vec q)}
```
on $u^{\alpha}_{i}$ representa el desplaçament de l'àtom $\alpha$ en la direcció $i$.

La segona derivada de l'energia de l'equació {eq}`Matriu_Dinàmica` correspon al canvi en la força que actua sobre l'àtom $\alpha$ en la direcció $j$ quan es desplaça l'àtom $\alpha'$ en la direcció $i$

\begin{equation}
\mathbf D^{\alpha',\alpha}_{i,j}=\frac{\partial}{\partial u^{\alpha*}_{i}}\vec F^{j}_{\alpha}(\vec q)
\end{equation}


```{math}
:label: matriu_dinàmica1
 \boxed{
 \mathbf D^{\alpha',\alpha}(\vec q)=\frac{1'}{\sqrt{M_\alpha' M_{\alpha}}}\sum_\vec l\vec\Phi^{\alpha',\alpha}\left(\vec R_\vec l-\vec R_{\vec l'}\right)e^{-i\vec q\cdot\vec R_\vec l}}
```


%Obtesses les posicions dels àtoms y classificats estos com primers, segons, etc. veïns, segons la distància al respectiu àtom de la ce\l.la $\vec 0$, procedim a calcular la contribució a la matriu dinàmica de cadascun dels àtoms, per la qual cosa necessitem conèixer el tensor de constants de for\c{c}a que correspon a la interacció de cada àtom amb el seu n-èssim veí.


Per calcular la matriu dinàmica {eq}`matriu_dinàmica1` emprant el model de constants de forces necessitem, per tant, constuir el tensor de constants de forces. Anem a suposar que les forces entre àtoms sols depenen del tipus d'elements químics que interaccionen i de la seua distància relativa.

Considerem un àtom $\alpha$ situat en la cel·la $\vec l$ , ($\alpha, \vec l$),a una certa distància, $|\vec R_{\alpha,\vec l}-\vec R_{\alpha',\vec{l'}}|$ de l'àtom  ($\alpha'$,$\vec l'$), i escollint el sistema de coordenades cartesianes de manera que la direcció del vector que uneix ambdós àtoms coincideix amb l'eix de les $x$, l'eix $y$ com la coordenada transversal en el plànol, $ti$, i $z$ la coordenada perpendicular al plànol $to$.

Podem escriure el tensor de forces d'este àtom, que segons el mòdul de la distància a l'àtom $(\alpha',\vec l$) classificarem com $n$-èssim veí, de la forma {cite:p}`wirtz04_phonon_disper_graph_revis`:


\begin{equation}
\mathbf \Phi_n^{\alpha',\alpha}=\begin{pmatrix}
\phi_{n,r}^{\alpha',\alpha}&\xi_n^{\alpha',\alpha} &0\\
-\xi_n^{\alpha',\alpha} & \phi_{n,ti}^{\alpha',\alpha} & 0 \\
0 & 0 & \phi_{n,to}^{\alpha',\alpha}
\end{pmatrix}
\label{eq:tensordeforces}
\end{equation}

L'estructura diagonal a blocs de la matriu reflexa el fet que estem supossant que  les vibracions interplanars y les de fora de pla, $to$, (en la direcció $z$) estan completament desacoblades.

Anem a suposar (simplificant encara més com són les interaccions entre àtoms) que un desplaçament longitudinal (radial, que estarà contés en el planol del cristall) o transversal en el planol, sols genera una força radial o transversal, respectivament, es a dir, $\xi_n^{\alpha,\alpha'}=0$ tal y com es realitza en {cite:p}`Balkanski_2000`.
%\missingfigure{Ací va imatge mostrant el cristall i les forces}


%\begin{figure}[h]
%\centering
%\includegraphics[width=40mm,height=40mm]{example-image-a}
%\caption{Caption}
%\end{figure}

Per tant, anem a considerar que el tensor de constants de forces  d'un àtom $\alpha$ classificat com $n$-èssim veí, situat en la direcció $\hat x$ respecte de l'átom  $\alpha'$, té la forma diagonal (notem que en la xarxa real que estem estudiant no té perque haver cap $n$-èssim vei en aquesta posició):

```{math}
:label: tensordeforcessimplificat
\mathbf\Phi_n^{\alpha'\alpha}=\begin{pmatrix}
\phi_{n,r}^{\alpha',\alpha}&0 &0\\
0& \phi_{n,ti}^{\alpha',\alpha} & 0 \\
0 & 0 & \phi_{n,to}^{\alpha',\alpha}
\end{pmatrix}
```


Per calcular el tensor de constants de forces de cadascun dels $i$ $n$-èssim veï real de l'àtom $\alpha'$, que formen un angle $\theta$ respecte de l'eix escollit com $x$, tenim que rotar la matriu de l'equació {eq}`tensordeforcessimplificat`:

\begin{equation}
 \mathbf\Phi_{n,i}^{\alpha'\alpha}(\theta)=\mathbf U^{-1}(\theta)\mathbf\Phi_n^{\alpha'\alpha}(0)\mathbf U(\theta)
\end{equation}


on $\mathbf U(\theta)$ ve donada per:
\begin{equation}
\mathbf U(\theta)=
\begin{pmatrix}
\cos(\theta)  & \sin(\theta) & 0 \\
-\sin(\theta) & \cos(\theta) & 0  \\
0             & 0            & 1
\end{pmatrix}
\end{equation}


Una vegada que sabem com tenim que construir el tensor de constants de forces, el càlcul de la matriu dinàmica és directe, ja que sols tenim que fer ús de l'equació {eq}`matriu_dinàmica1`, és a dir multipliquem el tensor de constants de forces associat a la interacció entre els àtoms ($\alpha', \vec l'=\vec 0$) i ($\alpha,\vec l$) per la fase, $e^{i \vec q\cdot \vec R_\vec l}$, on recordem, $\vec R_\vec l=l_1 \vec a_1+ l_2 \vec a_2$

En el cas considerat, com la base esta constituida per dos àtoms, la dimensió del tensor de constants de forces (i de la matriu dinàmica) és $3\cdot 2\times 3\cdot 2= 6\times 6$. Observem que esta matriu està escrita en termes de $4$ matris $3\times 3$:

\begin{equation}
 \label{eq:formamatriudinamica}
 \mathbf D=
 \begin{pmatrix}
  \mathbf D^{BB} && \mathbf D^{BN}\\
  \mathbf D^{NB} && \mathbf D^{NN}
 \end{pmatrix}
\end{equation}

recalquem que hem indexat el àtoms de la base pel tipus d'element al que pertanyen, pero en cas que els $2$ àtoms de la base pertanyeren al mateix tipus d'element químic la matriu dinàmica també seria $6x6$ (i podriem, per exemple, classificar-los com $A$ i $B$ per distingir-los).

Un punt important és que tenim que considerar (en les submatrius $D^{BB}$ i $D^{NN}$) les contribucions a la matriu dinàmica degudes a la interacció de l'àtom ($\alpha',0$) amb sí mateix. Ara bé, no necessitem escriure $\vec D_0^{BB}$ i $\vec D_0^{NN}$ explícitament, podem emprar el fet que si traslladem el conjunt d'àtoms en la mateixa direcció la força total és nul·la (si es desplaça el cristall com un tot no canvia l'energia potencial) i tal com ve indicat en {cite:p}`falkovsky08_symmet_const_phonon_disper_graph` a partir de l'equació {eq}`matriu_dinàmica1`, obtenim que (considerant fins a 3ers veïns):


\begin{equation}\begin{split}
\mathbf D_0^{BB}&=-\sum_{i=1}^3\mathbf \Phi_{1,i}^{BN}-\sum_{i=1}^6\mathbf \Phi_{2,i}^{BB}-\sum_{i=1}^3\mathbf \Phi_{3,i}^{BN} \\
\mathbf D_0^{NN}&=-\sum_{i=1}^3\mathbf \Phi_{1,i}^{NB}-\sum_{i=1}^6\mathbf \Phi_{2,i}^{NN}-\sum_{i=1}^3\mathbf \Phi_{3,i}^{NB}
\end{split}
\end{equation}


Es podrien tindre en compte altres simetries del cristall per determinar certes propietats del tensor de forces o de la seua transformada de Fourier, la matriu dinàmica, (certes relacions entre les components, tal com es fa en {cite:p}`falkovsky08_symmet_const_phonon_disper_graph`) però sols he emprat el fet que donat que l'energía es una funció quadràtica dels desplaçaments atòmics $u_{\alpha',\vec l'}, u_{\alpha,\vec l}$, la matriu de constants de forces té que ser simètrica,  i la seua transformada de Fourier, la matriu dinàmica té que ser una matriu hermítica, i per tant els seus valors propis, $\omega^2$ seran nombres reals.


Tenim per tant que la matriu dinàmica serà una matriu $6x6$ hermítica. Però, donat que les components en $z$ d'esta matriu es troben desacoblades, podem tractar aquestes vibracions de manera independent, de manera que en compte de treballar directament sobre una matriu $6x6$ de la forma

\begin{equation}
  \mathbf D(\vec q)=
  \begin{pmatrix}
   D_{xx}^{BB} & D_{xy}^{BB} & 0 &              & D_{xx}^{BN} & D_{xy}^{BN} & 0 \\
   D_{yx}^{BB} & D_{yy}^{BB} & 0 &              & D_{xx}^{BN} & D_{xy}^{BN} & 0 \\
       0       &     0       & D_{zz}^{BB} &    & 0  & 0  & D_{zz}^{BN}        \\
    &  &  &  &  &  \\
   D_{xx}^{NB} & D_{xy}^{NB} & 0 &              & D_{xx}^{NN} & D_{xy}^{NN} & 0 \\
   D_{yx}^{NB} & D_{yy}^{NB} & 0 &              & D_{xx}^{NN} & D_{xy}^{NN} & 0 \\
       0       &     0       & D_{zz}^{NB} &    & 0  & 0  & D_{zz}^{NN} 
  \end{pmatrix}
\end{equation}

treballarem sobre dues matrius 


 - Primer tractarem les vibracions transversal al pla del cristall mitjançant una matriu $2x2$:
 \begin{equation}
  \mathbf D_z(\vec q)=
  \begin{pmatrix}
    D_{zz}^{BB}  & D_{zz}^{BN}        \\
    D_{zz}^{NB}  & D_{zz}^{NN} 
  \end{pmatrix}
\end{equation}

 - Després tractarem les vibracions al pla del cristall emprant una matriu $4x4$:
 \begin{equation}
  \mathbf D_{xy}(\vec q)=
  \begin{pmatrix}
   D_{xx}^{BB} & D_{xy}^{BB} & D_{xx}^{BN} & D_{xy}^{BN}  \\
   D_{yx}^{BB} & D_{yy}^{BB} & D_{xx}^{BN} & D_{xy}^{BN}  \\
   D_{xx}^{NB} & D_{xy}^{NB} & D_{xx}^{NN} & D_{xy}^{NN} \\
   D_{yx}^{NB} & D_{yy}^{NB} & D_{xx}^{NN} & D_{xy}^{NN} \\
  \end{pmatrix}
\end{equation}
