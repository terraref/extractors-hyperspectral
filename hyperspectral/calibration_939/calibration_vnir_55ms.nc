�HDF

                    ���������3      ��������        `                               OHDR�"            ������������������������
               M                         title      �          Spectralon target with nominal visible reflectance = 0.95, as exposed to VNIR full image 1600 pixel and 268-298 lines on 20161021 ~13:15 local time in 2016_10_21_13_18_52_55ms. Spectralon is located in lines ~35-90 and samples (pixels) 600-1000..       created_by                   zenderT       nco_openmp_thread_number                                           l      history    L         Thu Jul 25 14:51:27 2019: ncks -O -d wavelength,0,938 calibration_vnir_55ms.nc mounted/calibration_vnir_55ms.nc
Mon Nov 14 08:42:30 2016: ncks -A -v xps_img_drk /Users/zender/data/terraref/clb/vnir_drk_avg_55ms.nc /Users/zender/data/terraref/clb/calibration_vnir_55ms.nc
Mon Nov 14 08:42:30 2016: ncks -O -v xps_img_wht /Users/zender/data/terraref/clb/vnir_wht_avg_55ms.nc /Users/zender/data/terraref/clb/calibration_vnir_55ms.nc
Mon Nov 14 08:26:27 2016: ncwa -O -a x,y /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_18_52_55ms/vnir_wht_cut_55ms.nc /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_18_52_55ms/vnir_wht_avg_55ms.nc
Mon Nov 14 08:26:27 2016: ncks -O -F -d x,600,1000 -d y,35,90 /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_18_52_55ms/vnir_wht_img_55ms.nc /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_18_52_55ms/vnir_wht_cut_55ms.nc
ncks -O --trr var_nm=xps_img_wht --trr_wxy=955,1600,269 --trr typ_in=NC_USHORT --trr typ_out=NC_USHORT --trr ntl_in=bil --trr ntl_out=bsq --trr ttl=Spectralon target with nominal visible reflectance = 0.95, as exposed to VNIR full image 1600 pixel and 268-298 lines on 20161021 ~13:15 local time in 2016_10_21_13_18_52_55ms. Spectralon is located in lines ~35-90 and samples (pixels) 600-1000. --trr_in=/Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_18_52_55ms/raw /Users/zender/terraref/computing-pipeline/scripts/hyperspectral/hyperspectral_dummy.nc /Users/zender/Downloads/VNIR_SpectralonRef_SinglePixel/2016_10_21_13_18_52_55ms/vnir_wht_img_55ms.nc �       NCO        `          netCDF Operators version 4.7.9 (Homepage = http://nco.sf.net, Code = http://github.com/nco/nco) �      history_of_appended_files          [         Mon Nov 14 08:42:30 2016: Appended file /Users/zender/data/terraref/clb/vnir_drk_avg_55ms.nc had following "history" attribute:
Mon Nov 14 08:28:45 2016: ncwa -O -a x,y /Users/zender/Downloads/VNIR-DarkRef/2016_10_19_04_17_07-55ms/vnir_drk_img_55ms.nc /Users/zender/Downloads/VNIR-DarkRef/2016_10_19_04_17_07-55ms/vnir_drk_avg_55ms.nc
ncks -O --trr var_nm=xps_img_drk --trr_wxy=955,1600,182 --trr typ_in=NC_USHORT --trr typ_out=NC_USHORT --trr ntl_in=bil --trr ntl_out=bsq --trr ttl=Dark counts as exposed to VNIR full image 1600 pixel and 182-218 lines on 20161019 ~3-4 AM local time in 2016_10_19_04_17_07-55ms. --trr_in=/Users/zender/Downloads/VNIR-DarkRef/2016_10_19_04_17_07-55ms/raw /Users/zender/terraref/computing-pipeline/scripts/hyperspectral/hyperspectral_dummy.nc /Users/zender/Downloads/VNIR-DarkRef/2016_10_19_04_17_07-55ms/vnir_drk_img_55ms.nc
���OHDR          �      �         !                       ���������           ������������������������0        CLASS                DIMENSION_SCALE     Q      n           A$      �                                    �wȐOCHK    �      N                              
wavelength�      ���OCHK`       NAME       @          This is a netCDF dimension but not a netCDF variable.       939 M�pOHDR          �      �                  
      ��      ��    �      v            ������������������������7     
   long_name                    Exposure counts          �       !       units                1  ��WOCHK    �      N                             xps_img_drk�      �L!OCHKH       meaning    (          Counts on scale from 0 to 2^16-1 = 655352       cell_methods       
          y, x: mean�s�OHDR          �      �                  
      ��      ��    �#      v            ������������������������7     
   long_name                    Exposure counts          �       !       units                1  �f��OCHK    #      z                             xps_img_wht�      9�~OCHKH       meaning    (          Counts on scale from 0 to 2^16-1 = 655352       cell_methods       
          y, x: mean�s�OCHK    %      V      P       DIMENSION_LIST                                             OqGCOL                        �                    �              �                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      OCHK     ������������������������J       _NCProperties      "          version=2,netcdf=4.6.2,hdf5=1.10.4Y��~                                                                            OCHK    m,      V      P       DIMENSION_LIST                                             ��h�OCHK�     t  REFERENCE_LIST       dataset                                       dimension                                                                      �              �              �{��                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           �P� j�P�l\�! �&b-'3:7f:R=�?�AJCyD�E�F�H�IrL�O�SoV�WdY ]�a#f�i�kn�pRtvv�x�y�z|�}�"���󈦉􉆌��Z��M����莦�}�F� ���j��;�-�d�o���K�`�C�x���Wʣ��Զ��׋����Z٦�ڂ���
�&�*�0�0�0�1�0�0�1�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�2�P�P�P�P�P�P�P�P�P�P�P�P�O�K���tʓ��#� ���������ؼ��5�E�K�L�N�O�N�N�N�M�N�N�N�L�J�J�I�B�>�>�;�+�$������ٯِ�_���ؖ�G��ט�Y�פ�Y���<�����k���ӝ�O�y��Π�K̤��X���Ȕ���|������#�f��B�C�������x��p�������+�㏦�G��|��0�?�؋�Аː������.�߇�b�I��8�Ň��݉���_�t���~�����K����"���N�%�㈄��j�L�_�K���~~�|�z�w�t�s�siu9w�x&y�x�w wvVu�t�sRsrp|m�j[hg�f�g�h2iFi�hch�g�f�eed@c�bb�a�`j`�_�^�]]�\D\\�[L[�Z�Y�XDXWW�VpV�UU�S�R�P�NXL�I�F�C�@>�;h9�7D6/5�456�7U9�:�;;x9�64�1�/�..�-�-�-B-V,f+�*�)�)o)X)�))*#+k,b-�-�-S-�,�+�*B*�)T)�('�$,"J�����yx��Lc$�+�4�0�0t������dD8AB,;|��]>b�h�����<��r2����_�-���qq��������		������K"���`