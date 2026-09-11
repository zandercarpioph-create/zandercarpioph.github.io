# Static preview server for the portfolio.
# Uses .NET HttpListener so the browser's idle preconnect sockets cannot
# stall the loop the way a raw TcpListener does.
param(
    [int]$Port = 8777,
    [string]$Root = "$PSScriptRoot\.."
)

$ErrorActionPreference = 'Stop'
$Root = (Resolve-Path $Root).Path

$types = @{
    '.html' = 'text/html; charset=utf-8'
    '.css'  = 'text/css; charset=utf-8'
    '.js'   = 'application/javascript; charset=utf-8'
    '.png'  = 'image/png'
    '.jpg'  = 'image/jpeg'
    '.jpeg' = 'image/jpeg'
    '.svg'  = 'image/svg+xml'
    '.pdf'  = 'application/pdf'
    '.json' = 'application/json'
    '.ico'  = 'image/x-icon'
}

$listener = New-Object System.Net.HttpListener
$listener.Prefixes.Add("http://localhost:$Port/")
$listener.Prefixes.Add("http://127.0.0.1:$Port/")
$listener.Start()
Write-Host "serving $Root on http://127.0.0.1:$Port/"

try {
    while ($listener.IsListening) {
        $ctx = $listener.GetContext()

        # One malformed or exotic request must never take the server down.
        try {
            $method = $ctx.Request.HttpMethod
            $isHead = ($method -eq 'HEAD')
            $path = [System.Uri]::UnescapeDataString($ctx.Request.Url.AbsolutePath)
            if ($path -eq '/') { $path = '/index.html' }

            $target = Join-Path $Root ($path.TrimStart('/') -replace '/', '\')
            $full = [System.IO.Path]::GetFullPath($target)

            # never serve anything outside the document root
            if (-not $full.StartsWith($Root, [StringComparison]::OrdinalIgnoreCase)) {
                $ctx.Response.StatusCode = 403
            }
            elseif (Test-Path -LiteralPath $full -PathType Leaf) {
                $bytes = [System.IO.File]::ReadAllBytes($full)
                $ext = [System.IO.Path]::GetExtension($full).ToLower()
                $ctx.Response.ContentType = if ($types.ContainsKey($ext)) { $types[$ext] } else { 'application/octet-stream' }
                $ctx.Response.ContentLength64 = $bytes.Length
                # HEAD carries the headers but no body
                if (-not $isHead) { $ctx.Response.OutputStream.Write($bytes, 0, $bytes.Length) }
                Write-Host "200 $method $path"
            }
            else {
                $body = [System.Text.Encoding]::UTF8.GetBytes("404 not found: $path")
                $ctx.Response.StatusCode = 404
                $ctx.Response.ContentType = 'text/plain; charset=utf-8'
                $ctx.Response.ContentLength64 = $body.Length
                if (-not $isHead) { $ctx.Response.OutputStream.Write($body, 0, $body.Length) }
                Write-Host "404 $method $path"
            }
        }
        catch {
            Write-Host "error: $($_.Exception.Message)"
        }
        finally {
            try { $ctx.Response.Close() } catch {}
        }
    }
}
finally {
    $listener.Stop()
    $listener.Close()
}
