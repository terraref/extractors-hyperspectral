�HDF

                    ���������3      ��������        `                               OHDR�"            ������������������������
               M                         title      �          Spectralon target with nominal visible reflectance = 0.95, as exposed to VNIR full image 1600 pixel and 268-298 lines on 20161021 ~13:15 local time in 2016_10_21_13_16_20_30ms. Spectralon is located in lines ~35-90 and samples (pixels) 600-1000..       created_by                   zenderT       nco_openmp_thread_number                                           l      history    L         Thu Jul 25 14:50:58 2019: ncks -O -d wavelength,0,938 calibration_vnir_30ms.nc mounted/calibration_vnir_30ms.nc
Mon Nov 14 08:42:30 2016: ncks -A -v xps_img_drk /Users/zender/data/terraref/clb/vnir_drk_avg_30ms.nc /Users/zender/data/terraref/clb/calibration_vnir_30ms.nc
Mon Nov 14 08:42:30 2016: ncks -O -v xps_img_wht /Users/zender/data/terraref/clb/vnir_wht_avg_30ms.nc /Users/zender/data/terraref/clb/calibration_vnir_30ms.nc
Mon Nov 14 08:26:09 2016: ncwa -O -a x,y /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_16_20_30ms/vnir_wht_cut_30ms.nc /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_16_20_30ms/vnir_wht_avg_30ms.nc
Mon Nov 14 08:26:09 2016: ncks -O -F -d x,600,1000 -d y,35,90 /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_16_20_30ms/vnir_wht_img_30ms.nc /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_16_20_30ms/vnir_wht_cut_30ms.nc
ncks -O --trr var_nm=xps_img_wht --trr_wxy=955,1600,296 --trr typ_in=NC_USHORT --trr typ_out=NC_USHORT --trr ntl_in=bil --trr ntl_out=bsq --trr ttl=Spectralon target with nominal visible reflectance = 0.95, as exposed to VNIR full image 1600 pixel and 268-298 lines on 20161021 ~13:15 local time in 2016_10_21_13_16_20_30ms. Spectralon is located in lines ~35-90 and samples (pixels) 600-1000. --trr_in=/Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_16_20_30ms/raw /Users/zender/terraref/computing-pipeline/scripts/hyperspectral/hyperspectral_dummy.nc /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_16_20_30ms/vnir_wht_img_30ms.nc �       NCO        `          netCDF Operators version 4.7.9 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco) �      history_of_appended_files          [         Mon Nov 14 08:42:30 2016: Appended file /Users/zender/data/terraref/clb/vnir_drk_avg_30ms.nc had following "history" attribute:
Mon Nov 14 08:26:55 2016: ncwa -O -a x,y /Users/zender/Downloads/VNIR-DarkRef/2016_10_19_03_00_27-30ms/vnir_drk_img_30ms.nc /Users/zender/Downloads/VNIR-DarkRef/2016_10_19_03_00_27-30ms/vnir_drk_avg_30ms.nc
ncks -O --trr var_nm=xps_img_drk --trr_wxy=955,1600,218 --trr typ_in=NC_USHORT --trr typ_out=NC_USHORT --trr ntl_in=bil --trr ntl_out=bsq --trr ttl=Dark counts as exposed to VNIR full image 1600 pixel and 182-218 lines on 20161019 ~3-4 AM local time in 2016_10_19_03_00_27-30ms. --trr_in=/Users/zender/Downloads/VNIR-DarkRef/2016_10_19_03_00_27-30ms/raw /Users/zender/terraref/computing-pipeline/scripts/hyperspectral/hyperspectral_dummy.nc /Users/zender/Downloads/VNIR-DarkRef/2016_10_19_03_00_27-30ms/vnir_drk_img_30ms.nc
F���OHDR          �      �         !                       ���������           ������������������������0        CLASS                DIMENSION_SCALE     Q      n           A$      �                                    �wȐOCHK    �      N                              
wavelength�      ���OCHK`       NAME       @          This is a netCDF dimension but not a netCDF variable.       939 M�pOHDR          �      �                  
      ��      ��    �      v            ������������������������7     
   long_name                    Exposure counts          �       !       units                1  ��WOCHK    �      N                             xps_img_drk�      �L!OCHKH       meaning    (          Counts on scale from 0 to 2^16-1 = 655352       cell_methods       
          y, x: mean�s�OHDR          �      �                  
      ��      ��    �#      v            ������������������������7     
   long_name                    Exposure counts          �       !       units                1  �f��OCHK    #      z                             xps_img_wht�      9�~OCHKH       meaning    (          Counts on scale from 0 to 2^16-1 = 655352       cell_methods       
          y, x: mean�s�OCHK    %      V      P       DIMENSION_LIST                                             OqGCOL                        �                    �              �                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      OCHK     ������������������������J       _NCProperties      "          version=2,netcdf=4.6.2,hdf5=1.10.4Y��~                                                                            OCHK    m,      V      P       DIMENSION_LIST                                             ��h�OCHK�     t  REFERENCE_LIST       dataset                                       dimension                                                                      �              �              �{��                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           Y,�?`������9��&�o!�"�#�$U%&�&�'K(�)�+�-J/0�0�2�5	8�9;w<�=�?%AbB�B�CODPEG�H�JtK�K LrM�N)PQ;Q�QQ	Q�N�K�IKvOW�[�^�_�a�e�hukk�i�h[i�kjo]s�wvzf}`Ȁ8�́��T�J�����אF��B�)�Ɣē��z����~����ߞ��N�����U�.��w��ۧ��٧S�o�e����Y�Z�٫٭C� ��_�̶����P�Y�_��� �B��{�Ñ���:��>���I���x�����#�}��r�«�{Ŕǧ�9��ɑ�d�(���#ɡ�3Ȗ���R�v�D���A��B�I�XҞ�n�)Ԕ����Ժ�Ԇ��ҩ�1���x�R�6�-�<���;��Ԗ���R֧��֩�֓����֌�7؀ؖ�Q�$��צף���!�^�_�O�<�=�����E؉���.�V�{يٌ�z�q�o�h�}نٮټ����������������ٌٛقٗٴ������������ٵ٩ٙ٤ٛ٦ٛ٢ٟپ���.�S�rڈڂ�wځ�{�v�o�t�sڐړڮڳ������ڻڭڒ�v�I�ڱ�K��'�r����E�a�g�n�vچڋڜڡڗڐړڑڇڐڥں������������������������������������ڰڥڠڤڦڶڸ������ڷڦڋ�l�_�F�<�6�&�*��'�#�+�D�T�c�uڑڏڟڷڶڶڿڽڿڳڣڑڇ�~�y�j�j�j�S�9�����ٯىفو٦ٻ������ٿق���7��ר��ר�H�����>�?�E�F�K�.�/�!������	�������������
���	���������������ٱ�{�=��������}�Q�� ĉōȷ�d̑���=�P�~����ϕ�E�.�qЅ�T�"��������̕��΅�>�f����� ��2�ѭ�m���φ�	�"�r�Y�g�Z��������[��v���&�ު����Ϥ����d���H�����e��ر������з÷�����h�Q�z�M�����.�P���ȸ�����߶,�\�>�����K�a�ƪj�&�ĩ�&����4�ږ~�|w<d<S�HE�FL�Q�X�`5j�s}����������u����ߔ��r�ⓡ�,�����א�{���C���������������C���"���h�.	~ }5|�{�zBz�y%y�xrxx�w.w�v�u�t�s�r
r�q�pGp�o�oWo�nnm�k�jihgteGcv`R]�YwV�SQ�NNiM�MHNAO>P�P�OhNbLlKvK�L�N�O�O�NcM=LCK�J�J�JhJXJ7JLJ�JK�K)L^L�L�L�L�L�L�L�LjLL�K�KmKSK+K�J�J�I%IH�FF�EEPD*C�A@J?V?C@EAB_B B�A&A�@5@�?j??^>H=�;=:�8;88�89j9w9A9�8�8 8�77�66�5l525�4�4	4�3 3�2s2;22�1�1M1�0J00�/S/9//�.9.�-�,�+�*m)�'[&�$#�!L !o���N �N  "�=��R���"�6�������,����:�X��T/�X�.�
9	��W#%7R�,#<	
n
N
�	D	�7���������������������������Lev����	|	�	d
�
�
�
�

4
�	s	'	���������		+	/	+	&	%	#				�����uX: