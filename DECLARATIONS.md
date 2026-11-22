# Declarations

## Acknowledgements and funding
The research presented in this paper was supported by the European Research
Council (ERC) under the European Union's Horizon Europe research and
innovation programme (grant agreement No. 101041711), by the Simons
Foundation as part of the Collaboration on the Mathematical and Scientific
Foundations of Deep Learning, by Heights Labs, by Convex Nexus Capital, by
the Israel Science Foundation (grant number 2258/19), and by the Israel
Science Foundation (ISF Grant 4101/25).

This work was done by multiple AI agents, Claude (Anthropic) and Codex (OpenAI), connected to Math Brain and working together inside the Brain framework, a shared, persistent memory in which the agents record claims, proofs, failed attempts and methods, retrieve one another's results, coordinate who works on what, and audit each other's work. For more on the Brain and Math Brain, see the Math Brain paper, the Math Brain repository (https://github.com/yspennstate/math-brain-public) and the Brain framework repository.

## Declaration of competing interest
The author holds equity in Heights Labs and in Convex Nexus Capital and
has a management role at Convex Nexus Capital. These entities are included
in the funding declaration above.

## Declaration of generative AI and AI-assisted technologies
The paper is part of Math Brain, a research project on AI systems proving mathematical results. In November 2026, Codex agents running OpenAI's GPT-6.1 proved alpha(S) <= I[f]/n for symmetric f and with it the bound Gamma <= 1.6 I, computed Gamma in closed form for threshold functions, and showed that Gamma/I tends to 1/ln 2 for majority. Claude Opus 5.5 (Anthropic) agents, with Claude Fable 5.2 as adviser, proved Theorems A and B, wrote and ran the computations of Section 6, searched the literature and prepared the manuscript. Other Claude agents re-derived the proofs by hand and checked them with independently written code. These checks are not peer review. The author directed the work and is responsible for the content of the paper and for decisions concerning its dissemination. The results rest on the work of others. The question is the Fourier entropy-influence conjecture of Friedgut and Kalai. The split of the entropy into a within-level part and a level part, and the first within-level bound, are due to O'Donnell, Wright and Zhou (Theorem 4 of their ICALP 2011 paper); the level bound used in the introduction is Lemma 2.1 of Shmalo, "Fourier entropy of threshold functions and their compositions", and the constant 2 + 1/ln 2 is from Shmalo, "The Fourier entropy-influence inequality for symmetric Boolean functions". The paper adds the exact within-level constant, the rate and the bound on Gamma.
