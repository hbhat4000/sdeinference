(TeX-add-style-hook
 "HMS2016"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("svmult" "graybox")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("footmisc" "bottom")))
   (TeX-run-style-hooks
    "latex2e"
    "svmult"
    "svmult10"
    "mathptmx"
    "helvet"
    "courier"
    "type1cm"
    "makeidx"
    "graphicx"
    "multicol"
    "footmisc"
    "amsmath"
    "amssymb"
    "latexsym")
   (LaTeX-add-labels
    "sec:1"
    "sec:2"
    "fig:pursuitdiagram"
    "eqn:pursuitODE"
    "eqn:pursuitSDE"
    "eqn:sde"
    "eqn:discretesde"
    "eqn:markovfactor"
    "subsec:2-1"
    "eqn:chapman"
    "eqn:phatalg"
    "subsec:2-2"
    "sec:3"
    "eqn:sdeapplication"
    "fig:lcosc"
    "fig:nbaspatial"
    "fig:inferredgamma")
   (LaTeX-add-bibliographies)))

