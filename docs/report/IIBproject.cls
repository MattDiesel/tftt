\NeedsTeXFormat{LaTeX2e}[1995/06/01]
\ProvidesClass{IIBproject}[2016/06/05 - Default margins 25mm. Title page changed. supervisor needed. Bugs to tl136@cam.ac.uk]
\newif\if@cuedwide
\@cuedwidefalse

\DeclareOption{a5paper}
    {\ClassError{IIBproject}{Option `a5paper' not supported}{}}
\DeclareOption{b5paper}
    {\ClassError{IIBproject}{Option `b5paper' not supported}{}}
\DeclareOption{twocolumn}%
    {\ClassError{IIBproject}{Option `twocolumn' not supported}{}}
\DeclareOption{titlepage}%
    {\ClassError{IIBproject}{Option `titlepage' not supported}{}}
\DeclareOption{wide}
    {\@cuedwidetrue\typeout{Wide option selected}}
\PassOptionsToClass{12pt,a4paper}{article}
\DeclareOption*{\PassOptionsToClass{\CurrentOption}{article}}
\ProcessOptions
\LoadClass{article}
\usepackage{textpos}
\usepackage{graphicx}
\usepackage{geometry}
\geometry{margin=25mm}
\if@cuedwide
\oddsidemargin=0pt \evensidemargin=0pt \topmargin=0pt \textwidth=421pt \textheight=612pt
\fi
\renewcommand{\abstractname}{Summary}
\def\supervisor#1{\gdef\@supervisor{#1}\gdef\thesupervisor{#1}}

\if@titlepage
  \renewcommand\maketitle{\begin{titlepage}%
  \let\footnotesize\small
  \let\footnoterule\relax
  \null\vfil
  \begin{textblock*}{120mm}(43mm,96mm)
  \Large \@title
  \end{textblock*}
  \@thanks
  \vfil\null
  \end{titlepage}%
  \setcounter{footnote}{0}%
  \let\thanks\relax\let\maketitle\relax
  \gdef\@thanks{}\gdef\@author{}\gdef\@title{}}
\else
\renewcommand\maketitle{\par
  \begingroup
    \renewcommand\thefootnote{\@fnsymbol\c@footnote}%
    \def\@makefnmark{\rlap{\@textsuperscript{\normalfont\@thefnmark}}}%
    \long\def\@makefntext##1{\parindent 1em\noindent
            \hb@xt@1.8em{%
                \hss\@textsuperscript{\normalfont\@thefnmark}}##1}%
    \if@twocolumn
      \ifnum \col@number=\@ne
        \@maketitle
      \else
        \twocolumn[\@maketitle]%
      \fi
    \else
      \newpage
      \global\@topnum\z@   % Prevents figures from going at top of page.
      \@maketitle
    \fi
    \thispagestyle{plain}\@thanks % tpl - remove page number 
  \endgroup
  \setcounter{footnote}{0}%
  \let\thanks\relax
  \let\maketitle\relax\let\@maketitle\relax
  \gdef\@thanks{}\gdef\@author{}\gdef\@title{}}
\def\@maketitle{%
  \newpage
  \null
  \includegraphics{Engineering.png}

  \begin{textblock*}{130mm}(33mm,43mm)
\LARGE \@title
\end{textblock*}
  \par

\begin{textblock*}{130mm}(33mm,122mm)
\begin{flushright}
Author Name: \@author{}\\
Supervisor: \@supervisor{}\\
\hspace{0cm}\\
Date: \today
\end{flushright}
\end{textblock*}

\begin{table}[b]
\begin{tabular}{l}
I hereby declare that, except where specifically indicated, the work\\ 
submitted herein is my own original work.\\
\\
Signed: \underline{\hspace{6cm}} Date: \underline{\hspace{3cm}}
\end{tabular}
\end{table}
\newpage}
\fi

\endinput

