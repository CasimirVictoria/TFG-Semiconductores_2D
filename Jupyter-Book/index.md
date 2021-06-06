Sobre este TFG
============================

## Objectius del treball proposat

  Els materials bidimensionals (2D) com el grafè són de gran interès tant per les seues propietats físiques exclusives com per les seues aplicacions potencials. L'estudi de la dinàmica de la xarxa cristallina (**fonons**) d'estos materials és un requisit previ per entendre la seua estabilitat estructural i propietats tèrmiques, així com les seues propietats de transport i òptiques.
  
  Este Treball de Fi de Grau consisteix en la computació dels modes vibracionals de materials semiconductors 2D y la seua correlació amb els observables rellevants per a la interpretació dels experiments de dispersió de la llum.

  L'estudiant treballarà la bibliografia pertinent per documentar-se sobre l'estat d'art teòric, i aprendrà una metodologia de càlcul basada en l'aproximació armònica a la dinàmica de una xarxa bidimensional. Seguidament, l'estudiant implementarà un model de constants de força per calcular els modes vibracionals d'un material 2D model.
Seguidament podran  determinar-se els modes vibracionals actius Raman.

Este treball està dissenyat especialment per a estudiants interessats en programació emprant llenguatges avançats tipus python, sempre orientats a aplicar l'aprés en les assignatures de Mecànica Quàntica o Física de l'Estat Sòlid.

## Metodologia

El càlculs s'han realitzar en python, emprant el software [sagemath](https://www.sagemath.org/index.html) i [jupyter](https://www.jupyter.org) com a entorn de treball. Els notebooks amb el càlculs es troben accessibles al repossitori públic de github [TFG-Semicunductores_2D](https://github.com/CasimirVictoria/TFG-Semiconductores_2D). També es pot consultar una versió en html del TFG, realitzada amb [jupyter{book}](https://jupyterbook.org/intro.html); en la versió online es mostra per defecte el codi emprat.

Jupyter/Jupyterlab ha demostrat ser un entorn fantàstic (i molt popular) per a programar en python (i altres llenguatges), però al no emprar fitxers de text pla (en este cas fitxers amb codi python), dificulta emprar eines de control de versions de codi com [git](https://git-scm.com) (entre altres coses). Per aquesta raó s'ha emprat també el plugin [jupytext](https://github.com/mwouts/jupytext) per a Jupyter, que permet tindre sincronitzat el notebook amb un fitxer en text pla amb el codi en python (com a fitxer python, markdown, ...), cosa que seu torn facilita enormement emprar el codi en altres entorns (editors de text, IDEs, ...), importar-lo, etc. 

