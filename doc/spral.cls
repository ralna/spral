\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{spral}[2014/03/06 SPRAL Documentation Class]

% Load base
\LoadClass{article}

% Library name definition
\RequirePackage{url}
\newcommand{\libraryname}{SPRAL}
\newcommand{\libraryurl}{\url{http://www.numerical.rl.ac.uk/spral}}

% Colors
\RequirePackage{pgf}
\definecolor{stfcblue}{RGB}{21,42,75}
\definecolor{stfcgrey}{RGB}{152,155,151}

% Required packages
\RequirePackage{fancyhdr}
\RequirePackage[top=2cm,left=2.5cm,right=2.5cm,bottom=2cm,headheight=12pt]{geometry}

% headers and footerrs
\fancypagestyle{firststyle}
{
   \fancyhf{}
   %\lfoot{\texttt{module \MakeLowercase{\libraryname\_\packagename}}}
   \lfoot{\libraryurl}
   \cfoot{\thepage}
   \rfoot{\versiondate}
   \renewcommand{\headrulewidth}{0pt}
}
\fancypagestyle{followstyle}
{
   \lhead{\bfseries \textcolor{stfcblue}{\libraryname} \texttt{\packagename}}
   %\lfoot{\texttt{module \MakeLowercase{\libraryname\_\packagename}}}
   \chead{Version \version}
   \lfoot{\libraryurl}
   \cfoot{\thepage}
   \rfoot{\versiondate}
}


% Beginning of document hooks
\AtBeginDocument{
   % Set headers
   \pagestyle{followstyle}
}

% Package information headers
\newcommand{\authors}{
}
\newcommand{\versionhistory}{
   \vspace{0.6cm}
   \noindent
   {\large \bfseries Major Version History}
}

% Support for specifying funders
\RequirePackage{tikz}
\newcommand{\funder}[2]{
   \hfill
   \begin{tikzpicture}
      \node (img) {\pgfimage[height=1cm]{#1}};
      \node[right,text width=3cm] at (img.east) {
         \footnotesize
         #2
      };
   \end{tikzpicture}
   \hfill
}

% New command definitions
\newcommand{\intentin}{\texttt{INTENT(IN)}}
\newcommand{\intentinout}{\texttt{INTENT(INOUT)}}
\newcommand{\intentout}{\texttt{INTENT(OUT)}}

% Redefine \maketitle
\makeatletter
\newcommand{\@pkglang}{}
\newcommand{\pkglang}[1]{
   \renewcommand{\@pkglang}{#1}
}
\def\@maketitle{
   \pgfimage[height=1cm]{stfc}
   \hfill
   {\Huge \bfseries \textcolor{stfcblue}{\libraryname}}
   \vspace{0.1cm}
   \textcolor{stfcgrey}{\hrule}
   \vspace{0.5mm}
   \hfill \textcolor{stfcgrey}{\Large \packagename\ v\version}
   \vspace{-0.3cm}

   \begin{center}
      {
         \LARGE \bfseries
         \texttt{\libraryname\_\packagename}\\
         \vspace{0.5em}
         \@title
      }

      \vspace*{0.5em}
      {\LARGE \@pkglang\ User Guide}
   \end{center}
   \begin{quote}
      \large
      \purpose
   \end{quote}

   \begin{flushright}
   %\noindent
   %{\bfseries Authors:}
   %\vspace{0.1cm}

   \noindent
   \@author
   \end{flushright}
   \vspace{-0.2cm}
   \textcolor{stfcgrey}{\hrule}
}
\makeatother

\newenvironment{funders}{
   \vspace{1cm}
   \textcolor{stfcgrey}{\hrule}
   \vspace{0.1cm}
}{
   \vspace{0.1cm}
   \textcolor{stfcgrey}{\hrule}
}

% Avoid hyphenation
\hyphenpenalty=10000
%\sloppy

\endinput
