from flask import Flask, render_template, request, redirect, url_for, flash, send_from_directory
import os
import time
import data_parser
import predictor
from google.cloud import storage

app = Flask(__name__)
app.secret_key = os.environ['SECRET_KEY']

# GCS
credentials = 'credentials.json'
storage_client = storage.Client.from_service_account_json(credentials)
bucket = storage_client.get_bucket('mzml-uploads')

def gcs_upload(file):
    try:
        timestamp = time.strftime("%Y%m%d-%H%M%S")
        blob = bucket.blob(timestamp)
        blob.upload_from_file(file.stream)
        return timestamp
    except Exception as err:
        print(err)

def gcs_download(blob):
    try:
        path = '/tmp/file.mzML'
        blob = bucket.blob(blob)
        blob.download_to_filename(path)
        return path
    except Exception as err:
        print(err)


@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        uploaded_file = request.files['file']
        if uploaded_file and uploaded_file.filename.split('.')[-1] == 'mzML':
            upload_name = gcs_upload(uploaded_file)
            download_path = gcs_download(upload_name)
            results = predictor.main(download_path)
            return redirect(url_for('results', r0=results[0], r1=results[1], r2=results[2], r3=results[3]))
        else:
            flash("File upload error. Please upload a .mzML file.")
    return render_template('index.html')


@app.route('/results')
def results():
    params = data_parser.run_parser('/tmp/file.mzML')
    results = (request.args.get('r0'), request.args.get('r1'), request.args.get('r2'), 
               request.args.get('r3'), params[0][0], params[1][1])
    return render_template('results.html', params=params, results=results)

@app.route('/download')
def download():
    return send_from_directory('/tmp', 'output.xml', as_attachment=True)


if __name__ == '__main__':
  app.run()
