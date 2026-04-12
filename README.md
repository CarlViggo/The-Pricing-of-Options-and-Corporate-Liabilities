# Analys av Early Exercise Premium: Europeiska vs. Amerikanska Optioner

Detta projekt undersöker prissättningen av amerikanska säljoptioner genom att lösa Black-Scholes PDE som en variationsolikhet med en finitadifferensmetod. Vi analyserar den "tidiga utnyttjandepremien" genom att numeriskt jämföra amerikanska och europeiska priser vid olika nivåer av volatilitet, ränta och löptid.

## Bakgrund

För en amerikansk säljoption $V(S,t)$ formuleras problemet som ett **Linjärt Komplementaritetsproblem (LCP)**. Optionens värde styrs av följande system:

$$
\begin{cases} 
\frac{\partial V}{\partial t} + \frac{1}{2}\sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + rS \frac{\partial V}{\partial S} - rV \leq 0 \\
V(S,t) \geq \max(K-S, 0) \\
\left( V - \max(K-S, 0) \right) \cdot \left( \frac{\partial V}{\partial t} + \mathcal{L}V \right) = 0
\end{cases}
$$

Genom att införa ett logiskt villkor i algoritmen kan den optimala utnyttjandegränsen $S^*(t)$ bestämmas. Denna gräns motsvarar punkten där det blir mer lönsamt att utnyttja optionen direkt än att behålla den, dvs. när räntevinsten på lösenpriset överstiger optionens tidsvärde.

Resultaten presenteras som relativa prisskillnader för att tydligt visa hur värdet av den amerikanska flexibiliteten varierar under olika marknadsförhållanden.

## Implementering
Projektet använder Crank–Nicolson-metoden för att uppnå stabila och noggranna beräkningar vid bestämning av den fria randen för det amerikanska utnyttjandet.