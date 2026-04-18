# Analys av Early Exercise Premium (EEP)

Detta projekt löser Black-Scholes PDE som en variationsolikhet för amerikanska säljoptioner:

$$
\begin{cases} 
\mathcal{L}_{BS}V \geq 0 \\
V \geq \max(K-S, 0) \\
(V - \text{payoff}) \cdot \mathcal{L}_{BS}V = 0
\end{cases}
$$

Implementeringen använder en Crank–Nicolson-metod där optionsvärdet vid varje tidsteg begränsas enligt:

$$V^{n+1} = \max\left( V_{pde}^{n+1}, K-S \right)$$

Metodbeskrivning finns här:

**[metodbeskrivning (PDF)](./Projektmetod.pdf)**