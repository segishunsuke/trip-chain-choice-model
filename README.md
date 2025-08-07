���̃��|�W�g���ł́C�g���b�v�`�F�C���S�̂̌��p���ő剻���闷�s�҂̈ӎv����𖾎��I�ɕ\�������C�ό����V�s�����f�����������߂̃v���O���������J���Ă��܂��B�g���b�v�`�F�C���̊ϑ��f�[�^��p�������f���̐����A���茋�ʂ𗘗p�����g���b�v�`�F�C���̗\�����s�����Ƃ��ł��܂��B

�v���O�����ɂ��Ă̖₢���킹�́C���L�̃��[���A�h���X�ɂ��肢���܂��D

<img src="./assets/images/mail.png" width="200px">


## Cython���C�u�����̃C���X�g�[��

���̃v���O������p����ɂ́CPython�������Cython���C�u�������K�v�ł��B

���g����Python����Cython���C�u�������C���X�g�[������Ă��Ȃ��ꍇ�́C�v�����v�g��ňȉ��̃R�}���h����͂��ăC���X�g�[�����s���ĉ������B

```
pip install Cython
```

### Windows���̒��ӓ_

Windows����Cython�v���O�������R���p�C������ɂ́CC/C++�̃R���p�C���iVisual Studio�j�̓������K�v�ł��B�R���p�C���͖����ŗ��p�\�ł��B

�ȉ��̎菇�ɏ]���āC�K�v�ȋ@�\���C���X�g�[�����Ă��������B

1. Visual Studio Community�̃C���X�g�[��

- [Microsoft�̌����T�C�g](https://visualstudio.microsoft.com/downloads/)�ɃA�N�Z�X���A�uVisual Studio Community�v���_�E�����[�h���܂��B
- �C���X�g�[�����N�����C�uPython�J���v��I�����ăC���X�g�[�����܂��B

2. Build Tools for Visual Studio�̃C���X�g�[��

- [�����y�[�W](https://visualstudio.microsoft.com/downloads/)�́uTools for Visual Studio�v�Z�N�V��������uBuild Tools for Visual Studio�v���_�E�����[�h���܂��B
- �C���X�g�[�����N�����C�uC++�ɂ��f�X�N�g�b�v�J���v��I�����ăC���X�g�[�����܂��B

## �v���O�����̃_�E�����[�h�ƃR���p�C��

[�ucodes�v�t�H���_](./codes)�t�H���_���̃t�@�C����S�ē���̃t�H���_�Ƀ_�E�����[�h���ĉ������B�t�@�C���̓��e�͈ȉ��̒ʂ�ł��B

- `trip_chain_simulator.pyx`, `trip_chain_simulator.pxd`�F�ό����V�s�����f�����������C�u�����̖{�̂ł��B
- `geneticr.pyx`, `geneticr.pxd`�F�����l�̈�`�I�A���S���Y���ɂ��֐��œK�����s�����C�u�����ł��B
- `mt19937ar.c`�F�����Z���k�c�C�X�^�ɂ�闐���������s��C�R�[�h�ł��B[�J���҂ɂ����J����Ă���R�[�h](https://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/MT2002/mt19937ar.html)��ҏW�������̂ł��B
- `setup.py`�F�R���p�C���ɗ��p����Python�t�@�C���ł��B
- `execute.py`�F`trip_chain_simulator`���C�u�����̗��p�Ⴊ�ڂ��Ă���Python�t�@�C���ł��B

�t�@�C����S�ă_�E�����[�h������A�ȉ��̃R�}���h�����s���ăv���O�������R���p�C�����ĉ������B

```
python setup.py build_ext --inplace
```

�R���p�C���͈�x���s����΁A����ȍ~�͍s���K�v�͂���܂���B




���̃v���O�������p����C�͏��W���̐ݒ���@���q�ׂ܂��D���̃v���O�����́C�L��`�P�f�ʂ����J���H�̕s�����v�Z�̊�b���ł���C
```math
\frac{dH}{dx} + \frac{1}{2g} \frac{d}{dx} \left( \frac{Q}{Bh} \right)^2 + \frac{n^2 Q^2}{B^2 h^{10/3}} = 0
```
��p���Ă��܂��D�����ŁC$`H`$�͐��ʕW��(m)��