PROJECT_ROOT=~/Documents/VUT/1MIT_zima/AVS/avs-proj02

rsync -azp ~/Documents/VUT/1MIT_zima/AVS/avs-proj02 barbora:~/
# ssh barbora 'cd ~/avs-proj02 && sbatch --wait evaluate.sl'
# ssh barbora 'ml intel && cd ~/avs-proj02/build_evaluate && make && ../scripts/compare.sh'

rsync -azp barbora:/home/ondrejmach/avs-proj02/build_vtune $PROJECT_ROOT/
rsync -azp barbora:/home/ondrejmach/avs-proj02/build_evaluate $PROJECT_ROOT/
rsync -azp --include='*.png' --exclude='*' barbora:/home/ondrejmach/avs-proj02/ $PROJECT_ROOT/
