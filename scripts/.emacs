(setq load-path  (cons (expand-file-name "~/.emacs.d/") load-path))

(require 'perltidy)
(add-hook 'cperl-mode-hook (lambda () (local-set-key (kbd "M-a") 'perltidy-dwim)))

;;(setq perl-tab-to-comment t)

(setq perl-indent-level 2)
(setq perl-continued-statement-offset 2)
(setq perl-continued-brace-offset -2)
(setq perl-tab-always-indent t)
(setq perl-brace-imaginary-offset 0)
(setq perl-label-offset -2)

(custom-set-variables
 '(cperl-indent-level 2)
 '(cperl-continued-statement-offset 2)
 '(cperl-tab-always-indent t)
 '(indent-tabs-mode nil)
)

(setq inhibit-startup-message t)
(setq-default transient-mark-mode t)

;; ========== Support Wheel Mouse Scrolling ==========
(mouse-wheel-mode t)

;; ========== Place Backup Files in Specific Directory ==========

;; Enable backup files.
(setq make-backup-files t)

;; Enable versioning with default values (keep five last versions, I think!)
(setq version-control t)

;; Automatically delete old versions.
(setq delete-old-versions t)

;; Save all backup file in this directory.
(setq backup-directory-alist (quote ((".*" . "~/.emacs_backups/"))))

;; show a menu only when running within X (save real estate when
;; in console)
(menu-bar-mode (if window-system 1 -1))

;; -----------------
;; Simple word count
;; -----------------

(defun word-count nil "Count words in buffer" (interactive)
  (shell-command-on-region (point-min) (point-max) "wc -w"))
