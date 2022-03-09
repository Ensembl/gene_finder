import logging

logging.basicConfig(filename='test_log.log', filemode='w', format='%(asctime)s,%(msecs)d %(name)s - %(levelname)s - %(message)s',datefmt='%H:%M:%S', level=logging.DEBUG)
def main():
    print('main')
    logging.info('train_model')
if __name__ == "__main__":
    try:
        logging.info("Start job Francesca ")
        
        main()
        
    except KeyboardInterrupt:
        
        
        print("Interrupted with CTRL-C, exiting...")
        sys.exit()