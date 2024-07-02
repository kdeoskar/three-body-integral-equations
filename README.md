# three-body-integral-equations

# Goal 
- To reproduce Figure 7 of https://arxiv.org/pdf/2010.09820. 

### Procedure
- Obtain $d_S^{(u,u)}$ numerically from Equation (34)
- Compute the residue $g$ of $\mathcal{M}_2$ at the bound-state pole using the following:
$$
    \begin{align*}
        s_b &= 4(m^2 - \kappa_{2k}^2) \text{        equation (17)} \\
        g &= 8\sqrt{2\pi \sqrt{s_b} \kappa_{2k}} \text{     equation (18)} \\
        \kappa &= 1/a \text{        (LO effective range expansion)}
    \end{align*}
$$
- Use $d_S^{(u,u)}$ and $g$ in Equation (24) to numerically calculate $\mathcal{M}_{\varphi b}$
$$ \lim_{\mathcal{K}_{df, 3} \rightarrow 0} \mathcal{M}_{\varphi b}(E) = g^2 \lim_{s_{2p}, s_{2k} \rightarrow s_{b}} d_S^{(u,u)} (p,k) $$

- The main goal of the code is to solve equation (34) numerically to obtain $d_S^{(u,u)}$.