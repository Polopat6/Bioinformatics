"""
notification_manager.py

Lightweight email + "text message" notifications for Monitor Mode (and,
by extension, any Advanced Mode background pipeline run) -- used to
alert a user to WARNINGS (e.g. a dropped folder had more FASTQ samples
than its metadata declared) and ERRORS (a folder was rejected, or a
launched pipeline run failed) without requiring them to keep a browser
tab open watching the Activity Feed.

--- Why plain smtplib (standard library), not a paid API ---
This app already avoids third-party network API dependencies wherever
a standard-library approach exists (see reference_manager.py's
urllib-based downloads, sra_manager.py's E-utilities calls). Every
organization running this portal already has (or can trivially set up)
an SMTP relay, so smtplib + STARTTLS covers "send an email" with zero
new package dependencies.

--- "Text message" notifications, via email-to-SMS gateways ---
Genuine SMS delivery normally requires a paid carrier/API (e.g.
Twilio), which isn't available as an installable dependency in this
environment. However, essentially every major US mobile carrier
operates a free email-to-SMS gateway: an email sent to
"<10-digit-number>@<carrier-gateway-domain>" is delivered as a text
message to that phone, with no account/API key needed on our end.
This module treats "send a text" as "send a very short email to this
gateway address" -- see CARRIER_SMS_GATEWAYS below. This is a
well-known, widely-documented technique, but it is NOT as reliable as
a real SMS API (some carriers rate-limit, delay, or occasionally
discontinue their gateway without notice) -- this caveat is surfaced
directly in the setup UI rather than silently promising guaranteed
delivery.

--- Configuration (deployment-level, shared across every monitor) ---
SMTP server/credentials are deployment-level environment variables,
NOT stored per-monitor, since these are almost always one shared
relay/account for the whole server:

    PIPELINE_SMTP_HOST            (required)
    PIPELINE_SMTP_PORT            (optional, default 587)
    PIPELINE_SMTP_USER            (optional, if auth required)
    PIPELINE_SMTP_PASSWORD        (optional, if auth required)
    PIPELINE_SMTP_FROM_ADDRESS    (required)
    PIPELINE_SMTP_USE_TLS         (optional, default "true")

Only the RECIPIENT address(es) -- email and/or phone+carrier -- are
part of a monitor's own per-monitor notification settings (see
monitor_manager.py's "notifications" config block). See
get_smtp_config_status() below for a fast ✅/⚠️ readiness check, shown
next to a monitor's "enable notifications" checkbox, mirroring
ontology_manager.check_gosemsim_available()'s pattern.

--- Never raises ---
A broken/misconfigured SMTP relay must NEVER crash a pipeline run or
the monitor daemon. Every public function here catches its own
exceptions and returns a (success, message) tuple (or a dict of them)
rather than propagating -- callers (monitor_manager.py's watch loop)
are expected to log the outcome to their own activity log, not treat a
failed notification as a reason to abort anything else.
"""
import os
import smtplib
import ssl
from email.mime.text import MIMEText
from email.utils import formatdate

# ---------------------------------------------------------------------------
# Common US carrier email-to-SMS gateways, for the "send as text" option.
# Not exhaustive, and not guaranteed to keep working forever (carriers can
# discontinue these without notice) -- offered as a free convenience, not a
# guaranteed-delivery channel. A user can always fall back to the email
# notification, which has no such dependency on carrier goodwill.
# ---------------------------------------------------------------------------
CARRIER_SMS_GATEWAYS = {
    "AT&T": "txt.att.net",
    "T-Mobile": "tmomail.net",
    "Verizon": "vtext.com",
    "Sprint": "messaging.sprintpcs.com",
    "Boost Mobile": "sms.myboostmobile.com",
    "Cricket Wireless": "sms.cricketwireless.net",
    "Google Fi": "msg.fi.google.com",
    "Metro by T-Mobile": "mymetropcs.com",
    "US Cellular": "email.uscc.net",
    "Xfinity Mobile": "vtext.com",
    "Visible": "vtext.com",
}


def build_sms_gateway_address(phone_number, carrier_label):
    """
    Build a "<digits>@<carrier-gateway-domain>" address for the given
    US phone number and a carrier label from CARRIER_SMS_GATEWAYS.
    Strips any non-digit characters from the phone number first, so
    "(503) 555-0100", "503-555-0100", and "5035550100" all work the
    same way; also tolerates an optional leading US country code "1".
    Returns None if the carrier isn't recognized or the number doesn't
    resolve to exactly 10 digits.
    """
    digits = "".join(c for c in str(phone_number) if c.isdigit())
    if len(digits) == 11 and digits.startswith("1"):
        digits = digits[1:]
    if len(digits) != 10:
        return None
    domain = CARRIER_SMS_GATEWAYS.get(carrier_label)
    if not domain:
        return None
    return f"{digits}@{domain}"


def _smtp_settings_from_env():
    return {
        "host": os.environ.get("PIPELINE_SMTP_HOST", ""),
        "port": int(os.environ.get("PIPELINE_SMTP_PORT", "587") or 587),
        "user": os.environ.get("PIPELINE_SMTP_USER", ""),
        "password": os.environ.get("PIPELINE_SMTP_PASSWORD", ""),
        "from_address": os.environ.get("PIPELINE_SMTP_FROM_ADDRESS", ""),
        "use_tls": os.environ.get("PIPELINE_SMTP_USE_TLS", "true").strip().lower() not in ("false", "0", "no"),
    }


def get_smtp_config_status():
    """
    Fast readiness check for whether notifications can actually be
    sent on this deployment, WITHOUT attempting a real SMTP
    connection -- mirrors ontology_manager.check_gosemsim_available()'s
    pattern: a cheap signal to show right next to an "enable
    notifications" checkbox, so a user knows immediately whether
    checking that box will do anything, before configuring a single
    monitor.
    Returns (configured: bool, message: str).
    """
    settings = _smtp_settings_from_env()
    if not settings["host"] or not settings["from_address"]:
        return False, (
            "Email/text notifications are not configured on this deployment "
            "-- an administrator needs to set the PIPELINE_SMTP_HOST and "
            "PIPELINE_SMTP_FROM_ADDRESS environment variables (see "
            "notification_manager.py's module docstring for the full list "
            "of supported variables)."
        )
    return True, f"✅ SMTP relay configured ({settings['host']}:{settings['port']}, from {settings['from_address']})."


def _send_email(to_address, subject, body):
    """
    Send a single plain-text email via the deployment's configured
    SMTP relay. Returns (success: bool, message: str). Never raises --
    every failure mode (missing config, connection error, auth error)
    is caught and returned as a clear message instead of propagating.
    """
    settings = _smtp_settings_from_env()
    if not settings["host"] or not settings["from_address"]:
        return False, "SMTP is not configured on this deployment (see get_smtp_config_status())."
    if not to_address:
        return False, "No recipient address was provided."
    msg = MIMEText(body)
    msg["Subject"] = subject
    msg["From"] = settings["from_address"]
    msg["To"] = to_address
    msg["Date"] = formatdate(localtime=True)
    try:
        with smtplib.SMTP(settings["host"], settings["port"], timeout=30) as server:
            if settings["use_tls"]:
                server.starttls(context=ssl.create_default_context())
            if settings["user"]:
                server.login(settings["user"], settings["password"])
            server.sendmail(settings["from_address"], [to_address], msg.as_string())
        return True, f"Sent to {to_address}."
    except (smtplib.SMTPException, OSError, TimeoutError) as e:
        return False, f"Failed to send to {to_address}: {e}"


def send_test_notification(email_address=None, sms_address=None):
    """
    Send a short "this works" test message to whichever address(es)
    are provided -- for a "Send test notification" button in the setup
    UI, so a user can confirm delivery BEFORE relying on it for an
    unattended, possibly multi-hour, monitor run.
    Returns {"email": (bool, str) or None, "sms": (bool, str) or None}.
    """
    results = {"email": None, "sms": None}
    subject = "Bioinformatics Portal -- Test Notification"
    body = (
        "This is a test notification from your Bioinformatics Portal's "
        "Monitor Mode notification settings. If you received this, "
        "notifications are configured correctly for this channel."
    )
    if email_address:
        results["email"] = _send_email(email_address, subject, body)
    if sms_address:
        # Kept deliberately short -- carrier SMS gateways often silently
        # truncate or drop long messages.
        results["sms"] = _send_email(sms_address, "", "Test alert: Bioinformatics Portal notifications are working.")
    return results


def _dispatch(notification_config, subject, body, sms_body=None):
    """
    Send subject/body to every enabled channel in notification_config
    (a plain dict shaped like monitor_manager.py's
    _effective_notification_config() output), silently no-op'ing any
    channel that isn't configured. Never raises.
    Returns a list of (channel, success, message) tuples, for the
    caller (monitor_manager.py) to fold into its own activity log.
    """
    if not notification_config or not notification_config.get("enabled"):
        return []
    outcomes = []
    email_address = notification_config.get("email_address")
    if email_address:
        success, message = _send_email(email_address, subject, body)
        outcomes.append(("email", success, message))
    sms_address = notification_config.get("sms_address")
    if sms_address:
        success, message = _send_email(sms_address, "", sms_body or body[:140])
        outcomes.append(("sms", success, message))
    return outcomes


def notify_warning(notification_config, monitor_id, folder_name, detail):
    """A non-blocking issue was detected (e.g. extra FASTQ files beyond what the metadata declares) -- the run still proceeded."""
    if not notification_config or not notification_config.get("notify_on_warning", True):
        return []
    subject = f"[Monitor: {monitor_id}] Warning -- {folder_name}"
    body = (
        f"Monitor '{monitor_id}' flagged a warning while processing the "
        f"dropped folder '{folder_name}':\n\n{detail}\n\n"
        "The pipeline run was still launched -- this is informational, "
        "not a blocking error."
    )
    return _dispatch(notification_config, subject, body, sms_body=f"[{monitor_id}] Warning on '{folder_name}': {detail[:100]}")


def notify_error(notification_config, monitor_id, folder_name, detail):
    """A dropped folder was rejected outright (missing FASTQ/metadata, or FASTQ/metadata sample counts didn't reconcile)."""
    if not notification_config or not notification_config.get("notify_on_error", True):
        return []
    subject = f"[Monitor: {monitor_id}] ERROR -- {folder_name}"
    body = (
        f"Monitor '{monitor_id}' could NOT process the dropped folder "
        f"'{folder_name}':\n\n{detail}\n\n"
        "This folder will not be retried automatically -- check its "
        "contents and, if needed, fix and rename/re-drop it so the "
        "watcher treats it as a new folder."
    )
    return _dispatch(notification_config, subject, body, sms_body=f"[{monitor_id}] ERROR on '{folder_name}': {detail[:100]}")


def notify_pipeline_error(notification_config, monitor_id, project_name, stage, message):
    """A launched background pipeline run failed partway through."""
    if not notification_config or not notification_config.get("notify_on_error", True):
        return []
    subject = f"[Monitor: {monitor_id}] Pipeline FAILED -- {project_name}"
    body = (
        f"The background pipeline run for project '{project_name}' "
        f"(launched by monitor '{monitor_id}') failed at stage "
        f"'{stage}':\n\n{message}\n\n"
        "Completed stages were preserved -- fixing the underlying issue "
        "and re-launching this project resumes from where it left off "
        "rather than starting over."
    )
    return _dispatch(notification_config, subject, body, sms_body=f"[{monitor_id}] Pipeline FAILED for '{project_name}' at '{stage}'.")


def notify_pipeline_success(notification_config, monitor_id, project_name):
    """A launched background pipeline run completed successfully."""
    if not notification_config or not notification_config.get("notify_on_success", False):
        return []
    subject = f"[Monitor: {monitor_id}] Pipeline complete -- {project_name}"
    body = (
        f"The background pipeline run for project '{project_name}' "
        f"(launched by monitor '{monitor_id}') completed successfully -- "
        "the gene counts matrix is ready. Open the project in the "
        "Differential Expression workspace to continue."
    )
    return _dispatch(notification_config, subject, body, sms_body=f"[{monitor_id}] Pipeline complete for '{project_name}'.")
