function time = tocf(mode) 
% TOCF ���Ԃ̕\���𐬌`�����X�g�b�v�E�I�b�`�^�C�}�̓ǂݍ���
%
%	TOCF(M)
%	TIC �� TOCF �֐��́A�o�ߎ��Ԃ𑪒肷�邽�߂Ɉꏏ�ɋ@�\���܂��B
%	TOC ���g�ł́ATIC ���g�p���Ă���̌o�ߎ��Ԃ𐮐��lM�Ŏw�肳�ꂽ�`���ŕ\�����܂��B
%	t = TOCF; �́A�o�ߓ����A���ԁA���A�b�̗v�f�������l�z��Ƃ��� t �Ɍo�ߎ��Ԃ�ۑ����܂��B
%
%	�\���`���ɂ��Ă͉��L�̕\���Q�l�ɂ��Ă��������B
%	�����Ƃ��ĕ\���`����^���Ȃ������ꍇ�́C1�`4�͈̔͂�0�̗v�f���Ȃ����̂������\�������`�����I������܂��B
%
%	---------------------------------------------------------------
%	| �����l |      �\���`��
%	---------------------------------------------------------------
%	|    1   |   �o�ߎ��Ԃ�11.123456�b�ł� 
%	|    2   |   �o�ߎ��Ԃ�12��11.123456�b�ł� 
%	|    3   |   �o�ߎ��Ԃ�1����12��11.123456�b�ł� 
%	|    4   |   �o�ߎ��Ԃ�1��1����12��11.123456�b�ł� 
%	|    5   |   Elapsed time is 12.123456 seconds. 
%	|    6   |   Elapsed time is 12 min 11.123456 sec. 
%	|    7   |   Elapsed time is 1 hr 12 min 11.123456 sec. 
%	|    8   |   Elapsed time is 1 day 1 hr 12 min 11.123456 sec. 
%	---------------------------------------------------------------
%
%	��
%	�X�g�b�v�E�H�b�`�̌o�ߎ��Ԃ𕪂ƕb�ŕ\�����܂��B
%	tocf(4)
%
%	�Q�l toc, tic, cputime.

% --
%	Title : TOCF()
%	Author : Sach1o : http://sach1o.blog80.fc2.com/
%	Created : 2007/12/12
% //--

% ���̓`�F�b�N�ƈ����Ȃ��̏ꍇ�̏���
if nargin==1 & [ mode<1 | mode>8 | ~isequal(mode, floor(mode))] 
	error('��1������1�`4�̐��̐����l�ł���K�v������܂��B');
end;

% ���A���A���A�b�ɕ�����
ntime=toc;
nday = floor(ntime/86400);
nhour = floor(mod(ntime,86400)/3600);
nmin = floor(mod(ntime,3600)/60);
nsec = mod(ntime,60);

% �ȉ��o�͏���
if nargout==1 
	time = [nday, nhour, nmin, nsec];
elseif nargout==0 
	if nargin==0 
		if nday>0 
			mode=4;
		elseif nhour>0 
			mode=3;
		elseif nmin>0 
			mode=2;
		else 
			mode=1;
		end;
	end;

	switch mode 
		case 1 
			strtime = sprintf( '%0.6f%s', ntime, '�b');
		case 2 
			nmin = nmin + 60*nhour + 1440*nday;
			strtime = sprintf( '%02d%s%04.6f%s', nmin, '��', nsec, '�b');
		case 3 
			nhour = nhour + 24*nday;
			strtime = sprintf( '%d%s%02d%s%04.6f%s', nhour, '����',nmin, '��', nsec, '�b');
		case 4 
			strtime = sprintf( '%d%s%d%s%02d%s%04.6f%s', nday, '��', nhour, '����',nmin, '��', nsec, '�b');
		case 5 
			strtime = sprintf( '%0.6f%s', ntime, ' seconds');
		case 6 
			nmin = nmin + 60*nhour + 1440*nday;
			strtime = sprintf( '%02d%s%04.6f%s', nmin, ' min ', nsec, ' sec');
		case 7 
			nhour = nhour + 24*nday;
			strtime = sprintf( '%d%s%02d%s%04.6f%s', nhour, ' hr ',nmin, ' min ', nsec, ' sec');
		case 8	 
			strtime = sprintf( '%d%s%d%s%02d%s%04.6f%s', nday, ' day ', nhour, ' hr ',nmin, ' min ', nsec, ' sec');
	end;

	if mode<5 
		disp(sprintf( '%s%s%s', '�o�ߎ��Ԃ�', strtime, '�ł�'));
	else 
		disp(sprintf( '%s%s%s', 'Elapsed time is ', strtime, '.'));
	end;
end;
