#!/bin/sh
#
# Completion for arts
#
# Author: Oliver Lemke  <olemke (at) core-dump.info>

shopt -s extglob progcomp

_artsexpand()
{
	[ "$cur" != "${cur%\\}" ] && cur="$cur"'\'

	# expand ~username type directory specifications
	if [[ "$cur" == \~*/* ]]; then
		eval cur=$cur
	elif [[ "$cur" == \~* ]]; then
		cur=${cur#\~}
		COMPREPLY=( $( compgen -P '~' -u $cur ) )
		return ${#COMPREPLY[@]}
	fi
}

_artsfiledir()
{
	local IFS=$'\t\n' xspec

	_artsexpand || return 0

	if [ "$1" = -d ]; then
		COMPREPLY=( ${COMPREPLY[@]} $( compgen -d -- $cur ) )
		return 0
	fi

	xspec=${1:+"!*.$1"}	# set only if glob passed in as $1
	COMPREPLY=( ${COMPREPLY[@]} $( compgen -f -X "$xspec" -- "$cur" )
		    $( compgen -d -- "$cur" ) )
}

_artsmethods()
{
    COMPREPLY=( ${COMPREPLY[@]} $( compgen -W '$( $1 -m all | tail -n +5 | \
		sed 's/^-//' | sed 's/^.-*.$//'  )' -- $cur ) )
}

_artsgroups()
{
	COMPREPLY=( ${COMPREPLY[@]} $( compgen -W '$( $1 -g | tail -n +5 | \
		sed 's/^-//' | sed 's/^.-*.$//'  )' -- $cur ) )
}

_artsvariables()
{
	COMPREPLY=( ${COMPREPLY[@]} $( compgen -W '$( $1 -w all | tail -n +5 | \
		sed 's/^-//' | sed 's/^.-*.$//'  )' \
		-- $cur ) )

}


_arts()
{
	local cur prev stopit

	COMPREPLY=()
	cur=${COMP_WORDS[COMP_CWORD]}
	prev=${COMP_WORDS[COMP_CWORD-1]}

	if [ $COMP_CWORD -gt 1 ]
	then
	    pprev=${COMP_WORDS[COMP_CWORD-2]}
	else
	    pprev=
	fi

	case "$pprev" in
	-@(w|-workspacevariables|i|-input|m|-methods|d|-describe))
	    return 0
	    ;;
	esac

    eval arts=$1
	case "$prev" in
	-@(w|-workspacevariables))
		_artsmethods $arts
		COMPREPLY=( ${COMPREPLY[@]} $( compgen -W 'all'  -- $cur ) )
		return 0
		;;
	-@(i|-input))
		_artsvariables $arts
		_artsgroups $arts
		return 0
		;;
	-@(I|-includepath|D|-datapath))
        COMPREPLY=( ${COMPREPLY[@]} $( compgen -d -- $cur ) )
        return 0
        ;;
	-@(m|-methods))
		_artsvariables $arts
		_artsgroups $arts
		COMPREPLY=( ${COMPREPLY[@]} $( compgen -W 'all' -- $cur ) )
		return 0
		;;
	-@(d|-describe))
		_artsvariables $arts
		_artsmethods $arts
		return 0
		;;
	-@(r|-reporting))
		COMPREPLY=( $( compgen -W '\
                 000 001 002 003 010 011 012 013 020 021 022 023 030 031 032 033 \
                 100 101 102 103 110 111 112 113 120 121 122 123 130 131 132 133 \
                 200 201 202 203 210 211 212 213 220 221 222 223 230 231 232 233 \
                 300 301 302 303 310 311 312 313 320 321 322 323 330 331 332 333 \
                  ' -- $cur ) )
		return 0
		;;
	-@(g|h|v|-groups|-help|-version))
		return 0
		;;
	*.arts)
		return 0
		;;
	*)
		_artsfiledir 'arts'
		;;
	esac


    if [ "$(echo $cur | sed 's/^\(.\).*$/\1/')" = "-" ]
	then
	# complete using basic options
	COMPREPLY=( $( compgen -W '--basename --datapath --describe --docdaemon \
                               --docserver --groups --help --input \
                               --includepath --methods --numthreads --plain \
                               --reporting --version --workspacevariables \
                               ' -- $cur ) )

	# this removes any options from the list of completions that have
	# already been specified somewhere on the command line.
	COMPREPLY=( $( echo "${COMP_WORDS[@]}" | \
	    (while read -d ' ' i; do
            [ "$i" == "" ] && continue
            # flatten array with spaces on either side,
            # otherwise we cannot grep on word boundaries of
            # first and last word
            COMPREPLY=" ${COMPREPLY[@]} "
            # remove word from list of completions
            COMPREPLY=( ${COMPREPLY/ ${i%% *} / } )
		done
		echo ${COMPREPLY[@]})
	    ) )
	fi

	_artsfiledir 'arts'
	return 0
}
complete -F _arts -o filenames arts

