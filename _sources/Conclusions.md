# Conclusions / *Conclusions*



- Observant els resultats que s'han obtés es pot concloure, que si bé el model de constants de forces és un model realment senzill i amb limitacions evidents, ja que podem observar que a altes freqüències el model és incapaç de reproduir la relació de dispersió calculada per primers principis emprant el mètode de *Density Functional Perturbation Theory*, este model sí és capaç de proporcionar una bona descripció de la relació de dispersió a baixes freqüències així com una descripció almenys qualitativa a altes freqüències (excepte per al mode LO on falla totalment) 

- Hem de tindre en compte que hem escollit el model més senzill possible, amb un tensor de constants de forces diagonal, i també que per realitzar l'ajust a les dades proporcionades sols s'han emprat $3$ punts dels més de $524$ disponibles (si bé són els punts crítics de la zona de Brillouin).
- Podem concloure, per tant, que este model de constants de forces és útils per obtindre un ajust ràpid si disposem de dades experimentals, i que amb quasi total seguretat pot millorar-se considerant més punts, per a la qual cosa no procediria obtenint les expressions analítiques als punts crítics com he fet ací, directament faria ús de numpy i scipy per realitzar l'optimització dels valors propis a les dades proporcionades).

- També és evident que quants més veïns considerem, millor ajust proporciona el model. Com podem comprobar en este treball l'ajust millora al quan considerem els quarts veïnss.

 - *Observing the obtained results we can conclude that the force-constants model is a really simple model with evident limitations: we can observe that at high frequencies the model is unable to reproduce the dispersion relation calculated by first principles using the *Density Functional Perturbation Theory* method, this model is capable of providing a good description of the  dispersion relation at low freqüències and a qualitative description at high freqüències (except for the  LO branch where totally fails)*

- We must keep in mind that we have chosen the simplest possible model, with a diagonal force-constant tensor, and also that to make the adjustment to the data provided  only $3$ points have been used out of the more than  $524$ available (although they are the critical points of the Brillouin area).

 - We can conclude, therefore, that this force-constants model is useful to obtain a quick adjustment if we have experimental data, and that with almost total certainty it can be improved considering more points, for which it I would not proceed calculating the analytical expressions at the critical points as I did here, I would directly make use of numpy and scipy to perform the optimization of the experimental eigenvalues.

- It is also evident that the more neighbors we consider, the better fit the model provides. As we can see in this work the adjustment improves when we consider the fourth neighbors.

