##################################
# Project Computer Graphics 2019 #
# Mano Marichal                  #
# Universiteit Antwerpen         #
##################################

Niet afgewerkte onderdelen:
- Stochastische L-systemen
- Buckyballs

Texture mapping:
Ik heb texture mapping voor bollen geimplementeerd. De ini files zijn textureXXX.ini

Uitbreidingen:
Ik heb 3 uitbreidingen:

* Moebiusbanden, de ini files heten Moebius.ini

* Vergelijkingen met 3 onbekenden:
Mijn graphics engine ondersteunt het tekenen van vergelijkingen met 3 onbekenden,
deze zijn van de vorm ax + by + cz + d = 0
De code loopt over alle X, Y en Z waardes binnen het interval [-range, range] met precision
Ik heb verschillende voorbeelden:
- euclidian_face001.ini: het X vlak bekenen vanuit 20, 20, 20
- euclidian_face002.ini: het x + y + z = 0 vlak bekenen vanuit 20, 20, 20
- euclidian_face003.ini: het X, Y en Z vlak bekenen vanuit 20, 20, 20
- euclidian_face004.ini: het XY, YZ en XZ vlak bekenen vanuit 20, 20, 20
- euclidian_face005.ini: random vlakken
- euclidian_face006.ini: Een hoop X = "constante" vlakken

* Projecteren van texturen als lichtbron
Bollen kunnen belicht worden met een textuur, wat dit doet is eigelijk per pixel de overeenkomstige texel van de textuur die de lichtbron uitzend te behandelen als zowel ambient, diffuus en speculair licht. Dit maakt het mogelijk om met schaduwen te werken, meerdere texturen op elkaar te projecteren en alleen bepaalde kleuren van de textuur te tekenen door je reflecties aan te passen.
De ini files zijn van de vorm texture_lightXXX.ini

