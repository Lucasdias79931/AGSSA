import os
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from email.mime.base import MIMEBase
from email import encoders


def send_email(recipient_email, zip_filepath, informacoes_path):
    smtp_server = os.getenv('MAIL_SERVER')
    smtp_sender = os.getenv('MAIL_USER')
    smtp_password = os.getenv('MAIL_SECRET')
    smtp_port = int(os.getenv('MAIL_ACCESS_PORT', 587))

    msg = MIMEMultipart()
    msg['From'] = smtp_sender
    msg['To'] = recipient_email
    msg['Subject'] = 'Resultado do Treinamento AGUA'

    body = "Seu processo foi concluído com sucesso. As informações do processamento são apresentadas a seguir e os resultados segue em anexo.\n\n"

    with open(informacoes_path, 'r') as f:
        body += f'{f.read()}\n'

    msg.attach(MIMEText(body, 'plain'))

    with open(zip_filepath, 'rb') as file:
        part = MIMEBase('application', 'zip')
        part.set_payload(file.read())
        encoders.encode_base64(part)
        part.add_header('Content-Disposition',
                        f'attachment; filename="{os.path.basename(zip_filepath)}"')
        msg.attach(part)

    try:
        server = smtplib.SMTP(smtp_server, smtp_port)
        server.ehlo()
        server.starttls()
        server.login(smtp_sender, smtp_password)
        server.sendmail(smtp_sender, recipient_email, msg.as_string())
        server.quit()

    except Exception as e:
        print(f'Erro ao enviar e-mail: {str(e)}')
