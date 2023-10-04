n_thread=8
python3 python/image_creation_improved.py 0 $(($1 / $n_thread)) $2 &
python3 python/image_creation_improved.py $(($1 / $n_thread)) $(($1 / $n_thread)) $2 &
python3 python/image_creation_improved.py $((2 * $1 / $n_thread)) $(($1 / $n_thread)) $2 &
python3 python/image_creation_improved.py $((3 * $1 / $n_thread)) $(($1 / $n_thread)) $2 &
python3 python/image_creation_improved.py $((4 * $1 / $n_thread)) $(($1 / $n_thread)) $2 &
python3 python/image_creation_improved.py $((5 * $1 / $n_thread)) $(($1 / $n_thread)) $2 &
python3 python/image_creation_improved.py $((6 * $1 / $n_thread)) $(($1 / $n_thread)) $2 &
python3 python/image_creation_improved.py $((7 * $1 / $n_thread)) $(($1 / $n_thread)) $2
