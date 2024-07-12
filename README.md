header-includes:
- \usepackage{ams}
- \usepackage{amsmath}
- \usepackage{amsfonts}

# three-body-integral-equations

# Goal 
- To reproduce Figure 7 of https://arxiv.org/pdf/2010.09820. 

### Procedure
1. Obtain $d_S^{(u,u)}$ numerically from Equation (34)
2. Compute the residue $g$ of $\mathcal{M}_2$ at the bound-state pole using the following:
- $$s_b = 4(m^2 - \kappa_{2k}^2) $$  equation (17)
- $$g = 8\sqrt{2\pi \sqrt{s_b} \kappa_{2k}} $$ equation (18)
- $$ \kappa = 1/a $$(LO effective range expansion)

3. Use $d_S^{(u,u)}$ and $g$ in Equation (24) to numerically calculate $\mathcal{M}_{\varphi b}$
$$ \lim_{\mathcal{K}_{df, 3} \rightarrow 0} \mathcal{M}_{\varphi b}(E) = g^2 \lim_{s_{2p}, s_{2k} \rightarrow s_{b}} d_S^{(u,u)} (p,k) $$

- The main goal of the code is to solve equation (34) numerically to obtain $d_S^{(u,u)}$.

### How do we actualy find the limit of $d_S^{(u,u)}$?
In order to evaluate $\lim_{\mathcal{K}_{df} \rightarrow 0} \mathcal{M}_{\varphi b}$ we need to find
$$ \lim_{s_{p}, s_{k} \rightarrow s_{b}} d_S^{(u,u)}(p, k) $$