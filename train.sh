CUDA_VISIBLE_DEVICES=0,1,2
chemprop_train --data_path ./data/GCN_cat_data_withN/GCN_data.csv --dataset_type multiclass --multiclass_num_classes 964 --save_dir ./data/GCN_cat_withN  --reaction --extra_metrics accuracy top3 --epochs 20
chemprop_train --data_path ./data/GCN_cat_data_withoutN/GCN_data.csv --dataset_type multiclass --multiclass_num_classes 963 --save_dir ./data/GCN_cat_withoutN   --reaction --extra_metrics accuracy top3 --epochs 20
chemprop_train --data_path ./data/GCN_solv_data_withN/GCN_data.csv --dataset_type multiclass --multiclass_num_classes 2973 --save_dir ./data/GCN_solv_withN  --reaction --extra_metrics accuracy top3 --epochs 20
chemprop_train --data_path ./data/GCN_solv_data_withoutN/GCN_data.csv --dataset_type multiclass --multiclass_num_classes 2972 --save_dir ./data/GCN_solv_withoutN  --reaction --extra_metrics accuracy top3 --epochs 20
chemprop_train --data_path ./data/GCN_reag0_data_withN/GCN_data.csv --dataset_type multiclass --multiclass_num_classes 4570 --save_dir ./data/GCN_reag0_withN  --reaction --extra_metrics accuracy top3 --epochs 20
chemprop_train --data_path ./data/GCN_reag0_data_withoutN/GCN_data.csv --dataset_type multiclass --multiclass_num_classes 4569 --save_dir ./data/GCN_reag0_withoutN  --reaction --extra_metrics accuracy top3 --epochs 20
