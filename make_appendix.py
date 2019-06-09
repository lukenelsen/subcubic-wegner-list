
from checking_realizations import *

from target_set_data import *

print_string = r'''\begin{longtable}{|l|l|}
\hline
\textit{Details} & \textit{Natural Core Subgraph}\\ \hline'''

for rc in NewTargetSet:  # c3a3, c3a4, c3a6, c4a4 are changed to cxa form in NewTargetSet in target_set_data.py
    NCS = NaturalCoreSubgraph(rc)
    if rc == "cxa33":
        NCS.coordinates[3] = (1,1)  # Centered and heightened a smidge.  Normal height is math.sin(math.radians(60)).
    print_string += r'''\hline

\hspace*{22mm}
\begin{minipage}{34mm}
\begin{itemize}
\item[\hspace*{5mm}Configuration:] \hypertarget{anchor'''+rc+r'}{'+DataDict[rc]['paper name']+r'''}

\item[\#Partitions:] '''+DataDict[rc]['partition count']+r'''

\item[\#Core Subgraphs:] '''+DataDict[rc]['core subgraph count']+r'''

\item[\#Realizations:] '''+DataDict[rc]['realization count']+r'''

\item[Runtime:] '''+DataDict[rc]['runtime']+r'''
\end{itemize}
\vspace*{2mm}
\end{minipage}

&
\begin{minipage}{115mm}

\centering

\begin{tikzpicture}[scale=1.56]
\begin{scope}[every node/.style={circle,draw},
                black/.style={shade,shading=axis,left color=black,right color=black,shading angle=90,inner sep=0pt,minimum size=8pt},]
                
'''
    # tikzpicture body.
    
    # Nodes.
    for v in NCS.graph.vertices():
        print_string += r'    \node ('+str(v)+r') at '+str(NCS.coordinates[v])+' [black] {};'+'\n'
    print_string += r'\end{scope}'+'\n\n'+r'\Large'+'\n'
    
    # Face length labels.
    FLL = NCS.facial_length_list
    for i in range(len(FLL)):
        if FLL[i] == 'x':
            print_string += r'    \node () at ('+str(i)+r'.5,0.2) {\*};'+'\n'
        else:
            print_string += r'    \node () at ('+str(i)+r'.5,0.2) {'+str(FLL[i])+'};'+'\n'
    if NCS.is_open:
        print_string += r'    \node () at ('+str(len(FLL)/2.0)+r',-0.3) {\*};'+'\n'
    else:
        print_string += r'    \node () at ('+str(len(FLL)/2.0)+r',-0.3) {'+str(NCS.central_face_length)+'};'+'\n'
    
    print_string += '\n'+r'\begin{scope}[every edge/.style={draw=black,very thick}]'+'\n'
    
    # Edges.
    for e in NCS.edges:
        print_string += r'    \path [-] ('+str(e[0])+r') edge node {} ('+str(e[1])+');'+'\n'
    # Special edge case:  If a configuration is closed and the central face length is exactly 1 greater than the number of edges in the spine, then the closing edge of the central face (the one not within the spine) will not show up because the spine is embedded as a horizontal line.  So we add another edge to curve it down.
    if not NCS.is_open and NCS.central_face_length == len(FLL)+1:
        out_str = str(45 - 5*(NCS.central_face_length-3))  # Adjust the angles to be thinner for longer chains.  (out_str will be preceded by a minus sign.)
        in_str = str(225 - 5*(NCS.central_face_length-3))
        print_string += r'    \path [-] (0) edge [out=-'+out_str+',in='+in_str+'] node {} ('+str(len(FLL))+');'+'\n'
    
    print_string += r'''\end{scope}
\end{tikzpicture}

\end{minipage}
\\ \hline'''

else:
    print_string += '\n\n'+r'\end{longtable}'






print print_string
