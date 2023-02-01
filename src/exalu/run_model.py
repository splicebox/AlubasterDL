import os
import time
import torch
import pickle

from .AluPrdEAD import AluPrdEAD
from .load_data import *
from .data_preprocess.curate_data import curate_data


def run_ead(strand, model_name='cnet',
            context_mode='none', single_side_pad_len=0,
            run_mode=0, work_dir=None, infer_set=None,
            infer_file=None, prd_file='prd_y.txt',
            model_wts_name=None,
            dataset_comment=None,
            model_dst_dir=None,
            dataset_save_file=None,
            test_chrs=['chr2', 'chr5']):
    since = time.time()
    print('test_chrs', test_chrs)
    # save path config
    if run_mode in [0, 1, 2, 4]:
        # run_mode 0, 2 will use their own dataset/model
        # run_mode 1 will only use 0's model, however, it will use its own dataset
        # run_mode 3 will use 2's dataset/model
        version_name = f'ead_v1.4t_{infer_set}'
        save_name_prefix = f'{version_name}_runmode_{run_mode}_padding_{context_mode}_{single_side_pad_len}'
        dataset_dst_dir = os.path.join(work_dir, f'dataset_pickles')
        if dataset_comment:
            # dataset comment is used for different dataset file, such as the snp part of code.
            dataset_save_file = os.path.join(dataset_dst_dir, f'{save_name_prefix}_dc_{dataset_comment}.pickle')
        else:
            dataset_save_file = os.path.join(dataset_dst_dir, f'{save_name_prefix}.pickle')
        if run_mode in [0, 2]:
            time_str = time.strftime('%Y%m%d_%H%M%S')
            model_dst_dir = os.path.join(work_dir, 'model_pts', f'{save_name_prefix}_{model_name}_{time_str}')

    print('dataset pickle file:', dataset_save_file)
    if run_mode == 0:
        # mode 0: only train
        neg_alu_fa = os.path.join(work_dir, 'neg_alu.fa')
        pos_alu_fa_file = os.path.join(
            work_dir, 'whole_pos_alu.fa')
        # load datasets
        if os.path.isfile(dataset_save_file):
            print('dataset pickle file exists')
            with open(dataset_save_file, 'rb') as fh:
                datasets = pickle.load(fh)
        else:
            print('dataset pickle file DOES NOT exist')
            curate_data(infer_set=infer_set, strand=strand, mode=context_mode, single_side_pad_len=single_side_pad_len, work_dir=work_dir)
            # datasets = load_data_ead_alu_chr_woinfer(alu_file=pos_alu_fa_file, neg_file=neg_alu_fa)
            datasets = load_data_ead_alu_chr_train_unique_duplicate(alu_file=pos_alu_fa_file, neg_file=neg_alu_fa, duplicate_times=10, test_chrs=test_chrs)
            os.makedirs(dataset_dst_dir, exist_ok='True')
            with open(dataset_save_file, 'wb') as fh:
                pickle.dump(datasets, fh)
        # train
        os.makedirs(model_dst_dir, exist_ok='True')
        alu_prd = AluPrdEAD(datasets, model_name=model_name, run_mode=run_mode)
        best_model_wts, model_records = alu_prd.train(run_mode, work_dir, infer_set, model_dst_dir)
        # save model
        torch.save(best_model_wts, os.path.join(model_dst_dir, 'best.pt'))
    if run_mode == 1:
        # mode 1: infer, use with mode 0, don't consider duplicates
        # the infer set is read by this mode self
        neg_alu_fa = os.path.join(work_dir, 'neg_alu.fa')
        if os.path.isfile(dataset_save_file):
            print('dataset pickle file exists')
            with open(dataset_save_file, 'rb') as fh:
                datasets = pickle.load(fh)
        else:
            print('dataset pickle file DOES NOT exist')
            datasets = load_data_ead_alu_chr_inferonly(strand, None, work_dir, infer_set)
            os.makedirs(dataset_dst_dir, exist_ok='True')
            with open(dataset_save_file, 'wb') as fh:
                pickle.dump(datasets, fh)
        alu_prd = AluPrdEAD(datasets, model_name=model_name, run_mode=run_mode)
        model_wts = torch.load(os.path.join(
            model_dst_dir, 'best.pt'))
        alu_prd.model.load_state_dict(model_wts)
        prd_y, y, id_line= alu_prd.evaluate('infer')
        prd_y = prd_y.tolist()
        y = y.tolist()
        assert(len(prd_y) == len(y) == len(id_line))
        with open(work_dir + '/' + prd_file, 'w') as write_fh:
            for i in range(len(id_line)):
                write_fh.write(f'{prd_y[i]}\t{y[i]}\t{id_line[i]}' + '\n')
    if run_mode == 2:
        # mode 2: train and infer every epoch
        neg_alu_fa = os.path.join(work_dir, 'neg_alu.fa')
        pos_alu_fa_file = os.path.join(
            work_dir, 'whole_pos_alu.fa')
        # load datasets
        if os.path.isfile(dataset_save_file):
            print('dataset pickle file exists')
            with open(dataset_save_file, 'rb') as fh:
                datasets = pickle.load(fh)
        else:
            print('dataset pickle file DOES NOT exist')
            curate_data(infer_set=infer_set, strand=strand, mode=context_mode, single_side_pad_len=single_side_pad_len, work_dir=work_dir)
            os.makedirs(dataset_dst_dir, exist_ok='True')
            datasets = load_data_ead_alu_chr_withinfer2(
                                            strand=strand,
                                            alu_file=pos_alu_fa_file,
                                            neg_file=neg_alu_fa,
                                            work_dir=work_dir,
                                            infer_set=infer_set,
                                            test_chrs=test_chrs,
                                            duplicate_times=10
                                            )
            with open(dataset_save_file, 'wb') as fh:
                pickle.dump(datasets, fh)
        # train
        os.makedirs(model_dst_dir, exist_ok='True')
        alu_prd = AluPrdEAD(datasets, model_name=model_name, run_mode=run_mode)
        best_model_wts, model_records = alu_prd.train(run_mode, work_dir, infer_set, model_dst_dir)
        # save model
        torch.save(best_model_wts, os.path.join(
            model_dst_dir, 'best.pt'))
    if run_mode == 3:
        # mode 3: infer, use with mode 2, consider duplicates, the infer set is from mode 2
        # load datasets
        if os.path.isfile(dataset_save_file):
            print('dataset pickle file exists')
            with open(dataset_save_file, 'rb') as fh:
                datasets = pickle.load(fh)
        alu_prd = AluPrdEAD(datasets, model_name=model_name, run_mode=run_mode)
        model_wts = torch.load(os.path.join(
            model_dst_dir, 'best.pt'))
        alu_prd.model.load_state_dict(model_wts)
        prd_y, y, id_line= alu_prd.evaluate('infer')
        prd_y = prd_y.tolist()
        y = y.tolist()
        assert(len(prd_y) == len(y) == len(id_line))
        with open(work_dir + '/' + prd_file, 'w') as write_fh:
            for i in range(len(id_line)):
                write_fh.write(f'{prd_y[i]}\t{y[i]}\t{id_line[i]}' + '\n')

    if run_mode == 4:
        # just the simplist infer
        # load datasets
        # if os.path.isfile(dataset_save_file):
        #     print('dataset pickle file exists')
        #     with open(dataset_save_file, 'rb') as fh:
        #         datasets = pickle.load(fh)
        # else:
        #     print('dataset pickle file DOES NOT exist')
        datasets = load_data_ead_simpleinfer(infer_file)
        # os.makedirs(dataset_dst_dir, exist_ok='True')
        # with open(dataset_save_file, 'wb') as fh:
            # pickle.dump(datasets, fh)
        alu_prd = AluPrdEAD(datasets, model_name=model_name, run_mode=run_mode)
        print(model_wts_name)
        model_wts = torch.load(model_wts_name)
        alu_prd.model.load_state_dict(model_wts)
        prd_y, y, id_line= alu_prd.evaluate('infer')
        prd_y = prd_y.tolist()
        y = y.tolist()
        assert(len(prd_y) == len(y) == len(id_line))
        with open(work_dir + '/' + prd_file, 'w') as write_fh:
            for i in range(len(id_line)):
                write_fh.write(f'{prd_y[i]}\t{y[i]}\t{id_line[i]}' + '\n')
    
    if run_mode == 5:
        # backward calculate gradients
        # load datasets
        if os.path.isfile(dataset_save_file):
            print('dataset pickle file exists')
            with open(dataset_save_file, 'rb') as fh:
                datasets = pickle.load(fh)
        alu_prd = AluPrdEAD(datasets=datasets, model_name=model_name, run_mode=run_mode)
        model_wts = torch.load(os.path.join(
            model_dst_dir, 'best.pt'))
        alu_prd.model.load_state_dict(model_wts)
        alu_prd.check_gradient()
    # time stamp
    time_elapsed = time.time() - since
    print('This run complete in {:.0f}m {:.0f}s'.format(
        time_elapsed // 60, time_elapsed % 60))
    
    return dataset_save_file, model_dst_dir
