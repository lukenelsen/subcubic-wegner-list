
from target_set_data import *

print_string = r'''\begin{table}[h!]
\begin{tabular}{|cc|}
\hline
\hspace*{5mm}&\\
&
\begin{minipage}{.95\textwidth}
\begin{multicols}{6}

'''

for rc in NewTargetSet:
    print_string += r'\hyperlink{anchor'+rc+'}{'+DataDict[rc]['paper name']+'}\n\n'

print_string += r'''\end{multicols}
\end{minipage}\\
&\\ \hline
\end{tabular}            
\caption{\label{tab:targetset}
The 76 reducible configurations used to prove Theorem~\ref{thm:prove}.
Each configuration is hyperlinked to a corresponding entry in the Appendix with more details.}
\end{table}'''

print print_string
