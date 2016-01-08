;;; dynare.el --- major mode to edit .mod files for dynare
;; Created: 2010 Sep 10
;; Version: 0.2

;; Copyright (C) 2010 Yannick Kalantzis
;;
;; This program is free software; you can redistribute it and/or modify
;; it under the terms of the GNU General Public License as published by
;; the Free Software Foundation, either version 3 of the License, or
;; (at your option) any later version.

;; This program is distributed in the hope that it will be useful,
;; but WITHOUT ANY WARRANTY; without even the implied warranty of
;; MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;; GNU General Public License for more details.
;;
;; You should have received a copy of the GNU General Public License
;; along with this program.  If not, see
;; <http://www.gnu.org/licenses/>.

;; Keywords: dynare

;; To comment/uncomment, use `ALT-;'. See `comment-dwim' for further details.

;;; Installation:
;;
;;   Put the this file as "dynare.el" somewhere on your load path, then
;;   add this to your .emacs or site-init.el file:
;;
;;   (require 'dynare)
;;   (autoload 'dynare-mode "dynare" "Enter dynare mode." t)
;;   (setq auto-mode-alist (cons '("\\.mod\\'" . dynare-mode) auto-mode-alist))

;;; Commentary:
;;
;;   Use C-c ; to comment/uncomment (or the default M-;).

;;; Changelog
;;
;;  2010-10-08 by Yannick Kalantzis v0.2
;;    Change syntax table: decimal numbers are words, mathematical
;;    operators are punctuation. `@' and `#' are symbols.
;;    Better font-lock for var, varexo, and parameters declaration,
;;    and for macrocommands.
;;    Add new keybinding for comment/uncomment: C-c ;
;;
;;  2010-09-07 by Yannick Kalantzis 
;;    Minor changes. Add `require'. Add builtin operators `;' and `='.
;;    Highlight lags and leads.
;;
;;  2010-09-06 by Yannick Kalantzis v0.1
;;    Created.  
;;    Reproduces Xah Lee's instructions.  See
;;    <http://xahlee.org/emacs/elisp_syntax_coloring.html>. Very basic
;;    syntax highlighting: comments, some keywords.

;;; TODO
;;    - indentation 
;;    - blocks templates "model/end", "initval/end", etc.
;;    - functions to insert main keywords

;;; Code:

;; function to comment/uncomment text
(defun dynare-comment-dwim (arg)
"Comment or uncomment current line or region in a smart way.
For detail, see `comment-dwim'."
   (interactive "*P")
   (require 'newcomment)
   (let ((deactivate-mark nil) (comment-start "//") (comment-end ""))
     (comment-dwim arg)))

;; define several class of keywords
(defvar dynare-keywords
  '("model" "initval" "endval" "end" "shocks" "periods" "values") 
  "dynare keywords.")

(defvar dynare-constructs-regexp
  "^[\s-]*\\(var\\|varexo\\|parameters\\)[a-zA-Z0-9,\s-@{}]*;"
  "regexp for dynare constructs.")

(defvar dynare-macrocommand-regexp
  "\\(@#define\\|@#include\\|@#endfor\\|@#for\\s-+[a-zA-Z0-9]+\\s-+in\\s-+[a-zA-Z0-9]+\\|@{[a-zA-Z0-9]+}\\)"
  "regexp for dynare macrocommand.")

(defvar dynare-functions
  '("simul" "stoch_simul" "steady" "check" "rplot" "dynatype" "dynasave" "resid")
  "dynare functions.")

;; create the regex string for each class of keywords
(defvar dynare-keywords-regexp (regexp-opt dynare-keywords 'words))
(defvar dynare-functions-regexp (regexp-opt dynare-functions 'words))

;; clear memory
(setq dynare-keywords nil)
(setq dynare-functions nil)

;; create the list for font-lock.
;; each class of keyword is given a particular face
(setq dynare-font-lock-keywords
  `(
    (,dynare-macrocommand-regexp . font-lock-builtin-face)
    (,dynare-functions-regexp . font-lock-function-name-face)
    (,dynare-keywords-regexp . font-lock-keyword-face)
    (,dynare-constructs-regexp . font-lock-keyword-face)
    (";\\|=" . font-lock-builtin-face)
    ("(\\(+\\|-\\)[1-9])" . font-lock-constant-face)
    ))

;; define the major mode
(define-derived-mode dynare-mode fundamental-mode
  "dynare mode"
  "dynare is a mode for editing mod files used by dynare."
  (setq mode-name "dynare mode")

  ;; modify the keymap
  (define-key dynare-mode-map [remap comment-dwim] 'dynare-comment-dwim)
  (define-key dynare-mode-map (kbd "C-c ;") 'dynare-comment-dwim)


  ;; decimal numbers should be treated as words
  (modify-syntax-entry ?\. "w" dynare-mode-syntax-table) 

  ;; mathematical operators are treated as punctuation
  (modify-syntax-entry ?+ "." dynare-mode-syntax-table)
  (modify-syntax-entry ?- "." dynare-mode-syntax-table)
  (modify-syntax-entry ?* "." dynare-mode-syntax-table)
  (modify-syntax-entry ?/ "." dynare-mode-syntax-table)
  (modify-syntax-entry ?^ "." dynare-mode-syntax-table)

  ;; symbols for the macrolanguage
  (modify-syntax-entry ?@ "_" dynare-mode-syntax-table)
  (modify-syntax-entry ?# "_" dynare-mode-syntax-table)

  ;; define C++ style comment  “/* ... */” and “// ...” 
  (modify-syntax-entry ?\/ ". 124b" dynare-mode-syntax-table) 
  ;; "/" is the 1st and 2nd char of /* and */ (a-style) and the 2nd char of //
  ;; (b-style)
  (modify-syntax-entry ?* ". 23" dynare-mode-syntax-table)
  ;; "*" is the 2nd and 1st char of /* and */ (a-style only)
  (modify-syntax-entry ?\n "> b" dynare-mode-syntax-table)
  ;; newline is the comment-end sequence of b-style comments

  ;; syntax highlighting
  (setq font-lock-defaults '((dynare-font-lock-keywords)))

  ;; clear memory
  (setq dynare-keywords-regexp nil)
  (setq dynare-functions-regexp nil)
  )

(provide 'dynare)
;;; dynare.el ends here
