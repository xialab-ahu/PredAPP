from flask import Flask,render_template,send_from_directory
from pathlib import Path
import time
from flask_mail import Mail,Message #邮件模块
from PredAPP import * 
import pandas as pd

from flask_wtf import FlaskForm
from wtforms import StringField,SubmitField,TextAreaField 
from flask_wtf.file import FileAllowed,FileRequired,FileField
from wtforms.validators import InputRequired,DataRequired,Length,ValidationError,Email

class SeqsForm(FlaskForm):
    sequences = TextAreaField(label='sequence',id='IDsequences',render_kw={'rows':'12','cols':'98%','class_':'layui-input-block'})
    email = StringField(label='mail',id='eMail',render_kw={'placeholder':'Enter your email','class_':'layui-input'})
    file = FileField(label='Sequence file',id='file',render_kw={"multiple data-min-file-count":'1'})
    submit = SubmitField('Submit',id='submitBtn',render_kw={'onclick':'check_sequence()','class_':'layui-btn'})

#flaskAPP设置部分
app = Flask(__name__, static_url_path='/PredAPP',static_folder='static')
app.secret_key = 'pRedAPP'
app.config['MAX_CONTENT_LENGTH'] = 5*1024*1024 
app.config['UPLOAD_PATH'] = Path(app.root_path).joinpath('/templates/uploads')
app.config['RESULT_PATH'] = Path(app.root_path).joinpath('/templates/pred_result')
Path(app.config['UPLOAD_PATH']).mkdir(exist_ok=True, parents=True)
Path(app.config['RESULT_PATH']).mkdir(exist_ok=True, parents=True)

#邮箱设置
app.config.update(
    MAIL_SERVER = 'smtp.qiye.aliyun.com', 
    MAIL_PORT = 465,
    MAIL_USERNAME = 'mail@xialab.info',
    MAIL_PASSWORD = 'your password',
    MAIL_DEFAULT_SENDER = ('PredAPP','mail@xialab.info')
)

mail = Mail(app)


@app.route('/') 
@app.route('/Index/')
@app.route('/PredAPP/')
def home():
    return render_template('home.html')


@app.route('/Results',methods=['GET','POST'])
@app.route('/PredAPP/Results',methods=['GET','POST'])
def Results():
    form = SeqsForm()
    seqs = form.sequences.data
    timestamp = time.strftime("%Y%m%d_%H%M%S",time.localtime())
    upload_seqs = Path(app.config['UPLOAD_PATH']).joinpath(str(timestamp)+'.txt')
    with open(upload_seqs,'wt') as f:
        f.write(seqs)


    resultfile = Model_pred(upload_seqs,timestamp,app.config['RESULT_PATH'])

    # 邮件发送部分
    if form.email.data:
        mail_addr = form.email.data
        mail_file_path = '{}/files/uploads/Mail/'.format(str(app.root_path))
        Path(mail_file_path).mkdir(exist_ok=True,parents=True)
        with open('{}/files/uploads/Mail/email.txt'.format(str(app.root_path)),'a') as f:
            f.write('{}\t{}\n'.format(mail_addr,timestamp))

        Subject = 'Predict Result - PredAPP' #发信主题
        Recipients = ['{}<{}>'.format(mail_addr,mail_addr)]#收件人
        Body = 'Your PredAPP analysis has been completed. Please refer to the attachment.' #邮件内容
        message = Message(subject=Subject,recipients=Recipients,body=Body)
        with app.open_resource(resultfile) as fp: #添加附件
            message.attach(filename='Pred_result.csv',content_type='text/plain', data=fp.read()) #filename为重命名附件的名字
        mail.send(message=message)

    #将结果csv文件以表格展示到网页
    df = pd.read_csv(resultfile)
    return render_template('result.html',tables=[df.to_html(header=True,table_id='table',classes='layui-table')],csvfile='/PredAPP/download/'+timestamp+'.csv')

@app.route('/download/<filename>',methods=['GET','POST']) #为生成的文件设置下载路径
@app.route('/PredAPP/download/<filename>',methods=['GET','POST']) #为生成的文件设置下载路径
def download(filename):
    down_path = app.config['RESULT_PATH']
    return send_from_directory(down_path,filename,as_attachment=True)

@app.route('/datasets/<filename>',methods=['GET','POST'])
def dataset(filename):
    return app.send_static_file("datasets/"+filename) 


@app.route('/Contact')
@app.route('/PredAPP/Contact')
def contact():
    return render_template('contact.html')

@app.route('/Help')
@app.route('/PredAPP/Help')
def help():
    return render_template('help.html')

@app.route('/RunClassifier')
@app.route('/PredAPP/RunClassifier')
def RunClassifier():
    form = SeqsForm() 
    return render_template('server.html',form = form) 

@app.route('/Download')
@app.route('/PredAPP/Download')
def download_page():
    return render_template('download.html')

#页面错误处理部分
@app.errorhandler(404)#定义404页面
def page_not_found(e):
    return render_template('errors/404.html'),404

@app.errorhandler(413)#定义413页面,文件大于5MB报错
def page_not_found(e):
    return render_template('errors/413.html'),413



if __name__ == '__main__':
    app.run()
